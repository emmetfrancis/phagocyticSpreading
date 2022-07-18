function outStruct = spreadCellCalcFEM_new(inStruct)
% this function calculates the next time step given the current shape of
% a spreading cell on a surface and spreading parameters.
% INPUTS: (these are contained as fields in the structure "inStruct")
% r0Adh: r0 parameter (range) in Lennard Jones potential (in nm)
% f0: scaling factor for Lennard Jones potential
% fPoly: maximum polymerization force (dimensionless - normalized by T0/d0)
% N: ligands per square micron
% dt: time step scaling factor (corresponds to desired Courant Number)
% mu: cytoplasmic viscosity in poise (10 Pa-s per poise)
% r: current (dimensionless) r coordinates of nodes (r/d0)
% z: current (dimensionless) z coordinates of nodes (r/d0)
% polyForceDecay: characteristic distance for decay of polymerization force
%   (given in microns)
% numIt: number of time steps
% 
% OUTPUTS: (in structure outStruct)
% rStored: cell of dimensionless r coordinates for each time point
% zStored: cell of dimensionless z coordinates for each time point
% v_rStored: cell of dimensionless r velocities for each time point
% v_zStored: cell of dimensionless z velocities for each time point
% pStored: cell of dimensionless pressure values for each time point
% timeVec: 1 by n vector of time values (in s)
% tensionProt: 2 by n vector - row 1 is dimensionless tension at each time point,
    % row 2 is protrusion strength (between 0 and 1) at each time point.

cellFig = gcf;
% Load inputs from structure
r0Adh = inStruct.r0Adh; % determines length scale of attraction force
f0 = inStruct.f0; % scaling factor
fPoly = inStruct.fPoly; % maximum protrusion stress
N = inStruct.N; % density of IgG in #/um^2
dtMult = inStruct.dt; % scaling factor for time step (in terms of Courant number)
dtRef = inStruct.dtRef; % reference value for max time
mu = inStruct.mu; % effective cytoplasmic viscosity
polyForceDecay = inStruct.polyForceDecay; % determines how much of the membrane the protrusive force acts over
numIt = inStruct.numIt; % total time = dt * numIt
d0 = inStruct.d0; % diameter of cell (defines dimensionless length)
r = inStruct.r;
z = inStruct.z;
volTol = inStruct.volTol; % percentage deviation in volume tolerated
discreteLogic = inStruct.discreteLogic;
NShift = inStruct.NShift; % r pos of first IgG molecule
% NShift = .5*sqrt(1/N);
% tempDecayFactor = 0;
muEff = mu*10; % possible re-scaling of mu
timeDecay = 8*muEff/2023;% 8 * muEff/2023; % for discrete adhesion case
% meshRef = inStruct.meshRef;
meshRefineTest = inStruct.meshRefineTest;
meshVal = inStruct.meshVal;
epsFactor = inStruct.epsFactor;
nodesPerEl = inStruct.nodesPerEl;
outStruct = struct; %  initialize output structure

% convert length parameters to dimensionless values
r0Adh = r0Adh/(d0*1e7);
% meshRef = r0Adh * meshRef;
meshVal = r0Adh * meshVal;
f0 = f0 * 1e-6 * 8.5^4; % new constant (8.5^4 relates to old one)
polyForceDecay = polyForceDecay/(d0*1e4);
N = N*(d0*1e4)^2;
% ligandDensity = N;
receptorDensity = 1000*8.5^2;
receptorFlux = (50*8.5^2);
contractStrength = 1;%.2;

% now set location of discrete adhesion sites for discrete simulations
if N > 0
    NGap = sqrt(1/N);
else
    NGap = 100;
end
if NShift == 0
    NLoc = 0;
else
    NLoc = [0,NShift/8.5];
end
% set locations of individual binding sites (could implement random dist
% here)
% inStruct.random = true;
while NLoc(end) < 5  % once beyond ~2, negligable contribution
    if inStruct.random
        curGap = random('Normal',NGap,NGap/2);
        while curGap < 0 || curGap > 2*NGap
            curGap = random('Normal',NGap,NGap/2);
        end
    else
        curGap = NGap;
    end
    nextLoc = NLoc(end) + curGap;
    NLoc = [NLoc, nextLoc]; % append locations one by one
end
isBound = false(1,length(NLoc));
bindLogic = false(1,length(NLoc));
adhDecayVal = 0;

% Define other parameters
maxAdhDist = 10*r0Adh; % how far from the surface to consider adhesion force
z0 = -r0Adh;
stationaryTime = 0;
areaList = nan(1,3);
aCoeff = [.05,.3]; % defines growth in polymerization stress
% aCoeff = [0, .4];
coeffCell = cell(1,2);
dx = 0.26; % determines smoothing of tension curve
transSA = 1.26;
transHeight = 0;
if meshRefineTest
    spacingVec = [.02,meshVal];
elseif fPoly > 0
    spacingVec = [.02,.05*r0Adh];
else
    spacingVec = [.02,.05*r0Adh];
end

% load in previous data (or starting data) and initialize storage cells
timeVec = inStruct.timeVec;
rStored = inStruct.rStored;
zStored = inStruct.zStored;
IENStored = inStruct.IENStored;
v_rStored = inStruct.v_rStored;
v_zStored = inStruct.v_zStored;
pStored = inStruct.pStored;
startEval = length(rStored)+1; % this is the first time step to solve for
rStored(startEval:startEval+numIt-1) = cell(1,numIt);
zStored(startEval:startEval+numIt-1) = cell(1,numIt);
v_rStored(startEval-1:startEval+numIt-2) = cell(1,numIt);
v_zStored(startEval-1:startEval+numIt-2) = cell(1,numIt);
pStored(startEval-1:startEval+numIt-2) = cell(1,numIt);
timeVec(startEval:startEval+numIt-1) = zeros(1,numIt);
contactAreaStored = zeros(1,startEval+numIt+1);
CaStored = zeros(1,startEval+numIt+1);
tensionProt = zeros(3,startEval+numIt+1);
errorLogic = zeros(1,startEval+numIt+1);
SA = zeros(1,startEval+numIt+1);
vol = zeros(1,startEval+numIt+1);
for i = 1:startEval-1
    curSurf = zStored{i} < eps;
    contactAreaStored(i) = pi*(max(rStored{i}(curSurf))*8.5)^2;
    [SA(i),vol(i)] = SAVolCalc(rStored{i},zStored{i},'trapz');
end
contactArea = contactAreaStored(startEval-1);

surfLogic = z <= eps;
z(surfLogic) = 0; % negative z points set to zero
startSurf = find(surfLogic,1,'first'); % eps threshold for numerical reasons
rStartSurf = interp1(z(1:startSurf),r(1:startSurf),0,'spline');
r(startSurf) = rStartSurf;
% for low densities, need to correct to have at least one segment of the
% cell adherent initially
if discreteLogic
    NLocVals = NLoc(NLoc < r(startSurf));
    isBound(NLoc < r(startSurf)) = true;
    if length(NLocVals) == 1 % correcting for any artifical folding back
        NLocVals = NLoc(1:2);
        isBound(1:2) = true;
        curIdx = startSurf-1;
        while r(curIdx) < NLocVals(2)
            r(curIdx) = NLocVals(2);
            curIdx = curIdx-1;
        end
    end
    r = [r(1:startSurf-1),NLocVals(end:-1:1)];
    z = [z(1:startSurf-1),zeros(1,length(NLocVals))];
end

% receptor model for discrete adhesion
N_Rtot = 1e6;
rho_Rcup = N_Rtot/(3*4*pi*(8.5/2)^2);
N_Rcup = N_Rtot * (contactArea/(3*4*pi*(8.5/2)^2));
rho_IgG = N/(8.5^2);
rho_Rbody = rho_Rcup;
N_Rbody = N_Rtot - N_Rcup;
K_a = 50;
bindArg = K_a + rho_Rcup + rho_IgG;
N_Rbound = (contactArea/2) * (bindArg - sqrt(bindArg^2 - 4*rho_Rcup*rho_IgG));
rho_Rfree = (N_Rcup - N_Rbound)/contactArea;
N_Rfree = N_Rcup - N_Rbound;
DBar = 1e-4 * sqrt(4*pi);
N_RStored = zeros(3,startEval+numIt+1);
N_RStored(:,startEval-1) = [N_Rbound; N_Rfree; N_Rbody];

% construct mesh for FEM calc
vMeshStruct = inStruct.vMeshStruct;
pMeshStruct = inStruct.pMeshStruct;
rStored{startEval-1} = vMeshStruct.xVec;
zStored{startEval-1} = vMeshStruct.yVec;
IENStored{startEval-1} = vMeshStruct.IEN;
if startEval == 2
    p0 = 4*ones(size(vMeshStruct.xVec)); % initial pressure should be ~4 (dimensionless) throughout (Law of Laplace)
else
    pInterp = scatteredInterpolant(rStored{startEval-2},zStored{startEval-2},pStored{startEval-2});
    p0 = pInterp(rStored{startEval-1},zStored{startEval-1});
end
r = vMeshStruct.xVec(vMeshStruct.boundIndices)';
z = vMeshStruct.yVec(vMeshStruct.boundIndices)';
z(z<eps) = 0; % correcting for any slight alt in contact due to interpolation
surfLogic = z<=eps;

% calculate s (arc length)
sParam = zeros(1,length(r));
for i = 2:length(r)
    sParam(i) = sParam(i-1) + sqrt((r(i)-r(i-1)).^2 + (z(i)-z(i-1)).^2);
end
SA0 = pi; % dimensionless initial surface area (4*pi*((d0/2)/d0)^2)
[SA(startEval-1),vol(startEval-1)] = SAVolCalc(r,z,'trapz'); % calculate current SA and volume
if SA(startEval-1) < (transSA-dx/2)*SA0
    T = .01 + .16*(SA(startEval-1)-SA0)/SA0;
elseif SA(startEval-1) < (transSA+dx/2)*SA0
    SAArg = (SA(startEval-1)-SA0)/SA0;
    T =  .7405*SAArg^2 - 0.0325*SAArg + 0.0225;
else
    T = .0516 + .54506*(SA(startEval-1)-transSA*SA0)/SA0;
end
T0 = .01; % resting cortical tension

%Initialize graph
plotContour(cellFig,vMeshStruct,z0,0)
k = startEval;
prevTime = timeVec(startEval-1);
while prevTime < (numIt-1)*dtRef
    
    %% PART 1: calculate cortex stresses
    [mCurv,pCurv,phiVals] = localCurvEval(r,z,false);
    pCurv(1) = mCurv(1); % from axisymm geometry
    pCurv(surfLogic) = 0; % curvature on the surface not considered
    mCurv(surfLogic) = 0;
    phiVals(surfLogic) = pi;
%     startAdh = find(z<5*r0Adh,1,'first');
%     mCurv(~surfLogic) = mean(mCurv(1:startAdh));
%     pCurv(~surfLogic) = mean(pCurv(1:startAdh));
    cStress_r = -(T/T0) * (mCurv + pCurv) .* sin(phiVals); % all stresses normalized by T0
    cStress_z = -(T/T0) * (mCurv + pCurv) .* cos(phiVals);
    
    %% PART 2: Calculate adhesion force
    adhStress_r = zeros(1,length(cStress_r));
    adhStress_z = zeros(1,length(cStress_z));
    adhLogic = ~surfLogic;
    rEdge = r(find(surfLogic,1,'first'));
    if f0 > 0 && N > 0
        if discreteLogic
            adhStress = adhForceDiscreteAxisymm_r4(r(~surfLogic),z(~surfLogic)-z0,maxAdhDist*10,r0Adh,f0,NLoc,rEdge,[]);
        else
            adhStress = adhForcePointDoubleTrapz_r4(r(~surfLogic),z(~surfLogic)-z0,maxAdhDist*10,r0Adh,f0,N,rEdge,N);
        end
    else
        adhStress = zeros(2,length(r(~surfLogic)));
    end
    adhStress = adhStress * exp(adhDecayVal);
    adhStress_r(adhLogic) = adhStress(1,:);
    adhStress_z(adhLogic) = adhStress(2,:);
    
    %% PART 3: Calculate polymerization force
    polyStress_r = zeros(size(cStress_r));
    polyStress_z = zeros(size(cStress_z));
    s0 = sParam(find(~surfLogic,1,'last')); % first point off surface -> max poly force
    polyLogic = ~surfLogic; % no poly force on surface
    SACur = SA(k-1);
    timePolyForce = 0; %overwritten if there is timePolyForce
    if fPoly > 0
        expCoeff = .5;
        if discreteLogic
            [timePolyForce,areaList,coeffCell] = timePolyCalc(SACur,contactArea,areaList,...
                expCoeff,aCoeff,coeffCell,0,transSA,dx);
        else
            [timePolyForce,areaList,coeffCell] = timePolyCalc(SACur,contactArea,areaList,...
                expCoeff,aCoeff,coeffCell,logical(f0)*N/8.5^2,transSA,0);
        end
        % alt simpler relationships for polymerization force over time
        % commented out below
%         timePolyForce = .1 + .9*(SACur-SA0)/(2*SA0);
%         timePolyForce = .02 + .98*(contactArea)/(250);
        if timePolyForce > 1
            timePolyForce = 1;
        end
        timePolyForce = timePolyForce / contractStrength;
        polyForce = fPoly * timePolyForce * exp((sParam(polyLogic)-s0)/polyForceDecay); % equal to fPoly at s0, less at earlier points
        polyForce = polyForce * exp(-stationaryTime/timeDecay); % decays the longer we go without a new contact (for discrete adhesion sites, for continuum case stationaryTime = 0)
        polyStress_r(polyLogic) = polyForce .* sin(phiVals(polyLogic));
        polyStress_z(polyLogic) = polyForce .* cos(phiVals(polyLogic));
    end
    tensionProt(:,k-1) = [T; timePolyForce; receptorDensity]; % save tension, relative polymerization stress, and receptor density at current time point
    
    %% PART 4: Solve FEM
    % formulate natural boundary condition (NBC)
    PList = vMeshStruct.P;
    IEN = vMeshStruct.IEN;
    boundIndices = vMeshStruct.boundIndices;
    hVecMesh = zeros(2*length(vMeshStruct.xVec),1);
    hVecMesh(2*boundIndices-1) = polyStress_r + adhStress_r + cStress_r;
    hVecMesh(2*boundIndices) = polyStress_z + adhStress_z + cStress_z;
    % now assemble hVec (using IEN for locating indices in linear system, cf. Hughes text)
    hVec = zeros(1,length(PList));
    hVec(1:2:end-1) = hVecMesh(2*IEN-1);
    hVec(2:2:end) = hVecMesh(2*IEN);
    
    vMeshStruct.hVec = hVec' .* logical(vMeshStruct.hVec); % only insert where NBC's are needed (where hVec is logic 1, see mesh fcn) 
    
    % now, solve linear system iteratively (Uzawa-type method)
    [v_r,v_z,KStruct,~,fBC_v] = FEMStokesVelocityCalc(vMeshStruct,p0,false,[],[]);
    [p,QStruct,~] = FEMStokesPressureCalc(pMeshStruct,v_r,v_z,p0,false,[],epsFactor);
    curError = max(abs(p-p0))/(max(abs(p))+min(abs(p))); % magnitude changes by (at most) curError*100%
    allError = curError;
    v_rFirst = v_r;
    v_zFirst = v_z;
    pFirst = p;
    numDiv = 0;
    dtAccelFactor = 1;
    errorLogic(k-1) = false;
    while abs(curError) > 1e-6
        p0 = p;
        [v_r,v_z] = FEMStokesVelocityCalc(vMeshStruct,p0,false,KStruct,fBC_v);
        p = FEMStokesPressureCalc(pMeshStruct,v_r,v_z,p0,false,QStruct,epsFactor);
        curError = max(abs(p-p0))/(max(abs(p))+min(abs(p)));
        allError = [allError,curError]; %#ok<AGROW>
        if length(allError) > 100 % not converging -> small time step and recalc
            v_r = v_rFirst;
            v_z = v_zFirst;
            p = pFirst;
            dtAccelFactor = 0.1;
            errorLogic(k-1) = 1;
            break
        end
            
        if ~(curError < allError(end-1))
            numDiv = numDiv + 1;
            if numDiv > 5 % not converging -> small time step and recalc
                v_r = v_rFirst;
                v_z = v_zFirst;
                p = pFirst;
                dtAccelFactor = 0.1;
                errorLogic(k-1) = 1;
                break
            end               
        end
    end
    v_rStored{k-1} = v_r;
    v_zStored{k-1} = v_z;
    pStored{k-1} = p;
    if any(isnan(v_r)) || any(isnan(v_z)) || any(isnan(p)) % throwing error
        error('some velocity or pressure is nan')
    end
        
    %% PART 5: Material point advection with surface contact threshold
    % isolate boundary velocities
    v_r = v_r(vMeshStruct.boundIndices);
    v_z = v_z(vMeshStruct.boundIndices);
    pSample = p(pMeshStruct.boundIndices);
    %convert from dimensionless velocities
    v_r = v_r' * (T0/muEff);
    v_z = v_z' * (T0/muEff);
    
    if fPoly == 0 % damping for adhesion only case
        startSurf = find(surfLogic,1,'first');
        idx = startSurf-5:startSurf-1;
        % compute stress balance at contact line, excluding fluid
        % velocities. If pressure is greater than the stress due to surface
        % tension, eqTerm < 0 and the adhesion stress is damped
        eqSample = pSample(idx)' - (T/T0) * (mCurv(idx)+pCurv(idx));
        eqTerm = trapz(sParam(idx), eqSample);
        if eqTerm < 0
            adhDecayVal = eqTerm/20; % dampening of adhesion stress
        else
            adhDecayVal = 0;
        end
    end
    
    % adaptive time step
    meshStackCur = cat(2,rStored{k-1},zStored{k-1},v_rStored{k-1},v_zStored{k-1});
    dtVal = dtCalc(meshStackCur,vMeshStruct.IEN,vMeshStruct.nodesPerEl); % dimensionless
    dt = dtMult * dtAccelFactor * dtVal * (d0/(T0/muEff));
    if dt < 1e-9 % indicates divergence
        errorLogic(k-1) = 1;
        timeVec(k) = timeVec(k-1) + dt;
        rStored{k} = [];
        zStored{k} = [];
        k = k+1;
        break
%         error('Simulation diverged (vanishing time step)')
    end
    % advect the points (simple advection first)
    r = r + v_r*dt/d0; % divided by d0 to convert back to dimensionless
    z = z + v_z*dt/d0;
    receptorDensity = receptorDensity + receptorFlux*dt;
%     curTime = timeVec(k-1);
%     if curTime > 220
%         MFIVal = 0;
%     else
%         MFIVal = interp1(CaTime,CaMFI,curTime);
%     end
%     timePolyForce = timePolyForce + MFIVal*(0.95/823.2)*dt;
    prevSurfLogic = surfLogic;
    sParam(1) = 0;
    for i = 2:length(r)
        sParam(i) = sParam(i-1) + sqrt((r(i)-r(i-1)).^2 + (z(i)-z(i-1)).^2);
    end
        
    % set points closer than a threshold to in contact with surface
    [SA(k),vol(k)] = SAVolCalc(r,z,'trapz'); % calculate current SA and volume
    if discreteLogic
%         bindLogic = (NLoc*8.5 - sqrt(contactArea/pi)) >= (1/sqrt(rho_Rbody));
%         diffusedLogic = (1/(DBar*sqrt(contactArea)*(rho_Rbody - rho_Rfree))) <= (tDiffusing + dt);
%         bindLogic = bindLogic | diffusedLogic;
        diffuseRange = sqrt(4*(DBar/sqrt(4*pi))*dt);
%         [r,z,~,errorFlag] = discreteThresh(r,z,NLoc,prevSurfLogic,T,SA(k),rho_Rbody*8.5^2,diffuseRange/8.5);
        [r,z,~,errorFlag] = discreteThresh(r,z,NLoc,prevSurfLogic,T,SA(k),inf,diffuseRange/8.5);
        if errorFlag
            break
        end
        surfLogic = z<=eps;
        startSurf = find(surfLogic,1,'first');
        isBound(NLoc < r(startSurf)) = true;
        curContactArea = pi*(8.5*r(startSurf))^2;
%         if curContactArea <= (contactArea+eps) && ~diffusedLogic
%             tDiffusing = tDiffusing + dt;
%         else
%             tDiffusing = 0;
%         end
        % advance receptor equation
        N_Rcup = N_Rcup + dt*(DBar*sqrt(curContactArea)*(rho_Rbody - rho_Rfree) +...
            ((curContactArea-contactArea)/dt)*rho_Rbody);
        rho_Rcup = N_Rcup/curContactArea;
        N_Rbody = N_Rtot - N_Rcup;
        rho_Rbody = (N_Rbody)/(3*4*pi*(8.5/2)^2 - curContactArea);
        bindingArg = K_a  + rho_Rcup + rho_IgG;
        N_Rbound = (curContactArea/2)*(bindingArg - sqrt(bindingArg^2 - 4*rho_Rcup*rho_IgG));
        N_Rfree = N_Rcup - N_Rbound;
        rho_Rfree = (N_Rfree)/curContactArea;
        N_RStored(:,k) = [N_Rbound; N_Rfree; N_Rbody];
        % compute stationary time (time since last binding to a new N)
        if abs(contactArea-curContactArea)<1e-12
            stationaryTime = stationaryTime + dt;
        else
            stationaryTime = 0;
        end
        if stationaryTime > timeDecay%2*(timeDecay/66)*muEff/2023
            break % then complete simulation (done moving)
        end
    else
        threshFraction = 0.0;
        newAttach = z <= -z0*threshFraction & ~prevSurfLogic; % points closer than threshold distance in z direction
        if any(newAttach)
            firstAttach = find(newAttach,1,'first');
            contactPt = z(firstAttach-1)/(z(firstAttach-1)-z(firstAttach)); % linear parameterization
            rContact = r(firstAttach-1) + contactPt * (r(firstAttach)-r(firstAttach-1));
            r(firstAttach) = rContact;
            z(newAttach) = 0;
        end
        surfLogic = z<=eps;
        startSurf = find(surfLogic,1,'first');
        curContactArea = pi*(8.5*r(startSurf))^2;
    end
    if curContactArea < contactArea-1e-12
        errorLogic(k-1) = nan;
        timeVec(k) = timeVec(k-1) + dt;
        rStored{k} = [];
        zStored{k} = [];
        k = k+1;
        break
%         error('retracted')
    end
    surfLogic(startSurf:end) = true; % everything after first contact is considered surface
    z(surfLogic) = 0;
    sParam = zeros(1,length(r));
    for i = 2:length(r)
        sParam(i) = sParam(i-1) + sqrt((r(i)-r(i-1)).^2 + (z(i)-z(i-1)).^2);
    end
    maxGap = max(diff(sParam));
    span = round(.2/(2*maxGap))*2 + 1; % span for smoothing
    rSmooth = smooth(sParam(1:startSurf-4),r(1:startSurf-4),span,'loess');
    zSmooth = smooth(sParam(1:startSurf-4),z(1:startSurf-4),span,'loess');
    negLogic = zSmooth <= eps; % in case any smoothing advects points to "contact"
    rSmooth(negLogic) = r(negLogic);
    zSmooth(negLogic) = z(negLogic);
    r(1:startSurf-4) = rSmooth;
    z(1:startSurf-4) = zSmooth;
    r(1) = 0;
    [r,z] = volDilation(r,z,pi/6,volTol,false); % corrects for volume deviations by normal displ of boundary nodes
    
    %% PART 8: Respace nodes and then smooth contour and correct volume
    [r,z] = respaceNodes(r,z,transHeight,spacingVec); % respace function
    if pi*(8.5*r(find(z<=0,1,'first')))^2 < curContactArea % if respacing causes retraction, error (should not happen)
        error('retracted')
        %             r(find(z<=0,1,'first')) = sqrt(curContactArea/pi)/8.5;
    end
    % recalculate arc length
    sParam = zeros(1,length(r));
    for i = 2:length(r)
        sParam(i) = sParam(i-1) + sqrt((r(i)-r(i-1)).^2 + (z(i)-z(i-1)).^2);
    end
    
    surfLogic = z <= eps; %logical vector of surface nodes
    startSurf = find(surfLogic,1,'first');
    contactRad = r(startSurf);
    contactArea = pi*(d0*1e4)^2*contactRad^2;
    
    %% PART 9: Saving data, update SA, vol, and tension, then display contour
    [vMeshStruct,pMeshStruct] = meshFromBound(r,z,nodesPerEl);
    rStored{k} = vMeshStruct.xVec;
    zStored{k} = vMeshStruct.yVec;
    IENStored{k} = vMeshStruct.IEN;
    pInterp = scatteredInterpolant(rStored{k-1},zStored{k-1},pStored{k-1});
    p0 = pInterp(rStored{k},zStored{k});
    r = vMeshStruct.xVec(vMeshStruct.boundIndices)';
    z = vMeshStruct.yVec(vMeshStruct.boundIndices)';
    z(z<eps) = 0;
    surfLogic = z<=eps;

    % SA-tension relationship determined in Herant et al 2005
    [SA(k),vol(k)] = SAVolCalc(r,z,'trapz'); % calculate current SA and volume
    contactAreaStored(k) = contactArea;
%     T = -.24 + (.25*contractStrength)*exp((SA(k)-SA0)/SA0);
    if SA(k) < (transSA-dx/2)*SA0
        %     T = .025;
        T = .01 + .16*(SA(k)-SA0)/SA0;
    elseif SA(k) < (transSA+dx/2)*SA0
        %     T = (.1/(2*dx))*(SA(startEval-1)/SA0-(transSA-dx/2))^2 + .025;
        SAArg = (SA(k)-SA0)/SA0;
        T =  .7405*SAArg^2 - 0.0325*SAArg + 0.0225;
    else
        %     T = .025 + 0.1 * (SA(startEval-1)/SA0-transSA);
        T = .0516 + .54506*(SA(k)-transSA*SA0)/SA0;
    end
    timeVec(k) = timeVec(k-1) + dt;
    plotContour(cellFig,vMeshStruct,z0,timeVec(k)) % update graph
    prevTime = timeVec(k);
    k = k+1;
    sParam = zeros(1,length(r));
    for i = 2:length(r)
        sParam(i) = sParam(i-1) + sqrt((r(i)-r(i-1)).^2 + (z(i)-z(i-1)).^2);
    end
end

% save to outStruct
if length(zStored) < k-1
    k = length(zStored) + 1;
end
outStruct.rStored = rStored(1:k-1);
outStruct.zStored = zStored(1:k-1);
outStruct.IENStored = IENStored(1:k-1);
outStruct.vol = vol(1:k-1);
outStruct.SA = SA(1:k-1);
outStruct.timeVec = timeVec(1:k-1);
if discreteLogic
    outStruct.NLoc = NLoc;
    outStruct.N_RStored = N_RStored(:,1:k-1);
end
outStruct.pStored = pStored(1:k-2);
outStruct.v_rStored = v_rStored(1:k-2);
outStruct.v_zStored = v_zStored(1:k-2);
outStruct.tensionProt = vertcat(tensionProt(:,1:k-2),errorLogic(:,1:k-2));
end


%% Define other functions
function [] = plotContour(cellFig,meshStruct,z0,timeVal)
xVec = meshStruct.xVec;
yVec = meshStruct.yVec;
nodesPerEl = meshStruct.nodesPerEl;
figure(cellFig)
cla
hold on
title(sprintf('t = %.2f s',timeVal))
switch nodesPerEl
    case 4
        numCols = meshStruct.numCols;
        numRows = meshStruct.numRows;
        xMesh = reshape(xVec,[numCols,numRows])';
        yMesh = reshape(yVec,[numCols,numRows])';
        surf(xMesh,yMesh,zeros(size(xMesh)),'FaceColor','none','EdgeColor','k')
        surf(-xMesh,yMesh,zeros(size(xMesh)),'FaceColor','none','EdgeColor','k')
        view([0 90])
    case 3
        IEN = reshape(meshStruct.IEN,[3,length(meshStruct.IEN)/3])';
        triplot(IEN,xVec,yVec,'Color','k')
        triplot(IEN,-xVec,yVec,'Color','k')
end
%         plot(r, z, 'r*')
%         plot(-r, z, 'r*')
x = [-1.5, -1.5, 1.5, 1.5];
y = [0, -z0, -z0, 0];
patch('XData',x,'YData',y','FaceColor','k')
xlim([-1.5 1.5])
ylim([z0 1.5])
daspect([1 1 1])
drawnow
end

function [timePolyForce, areaList, coeffCell] = timePolyCalc(SACur,contactArea,areaList,n,aCoeff,coeffCell,N,transSA,dx) 
    xVecData = [98.1575363494768,119.818133608253;
                98.3382826750535,119.960330704836;
                98.6084670863954,120.331727399340;
                99.3756925868717,121.088347532733;
                100.521205134533,122.500829954456;
                102.011184064969,124.342244628163]; % from testing, known contact areas associated with certain deformations
               % col 1 is CA for 26% SA deformation, col 2 is CA for 39%
    NVals = [0,100,300,1000,3000,10000];
    xVec = [0,0,0];
    xVec(2) = interp1(NVals,xVecData(:,1),N);
    xVec(3) = interp1(NVals,xVecData(:,2),N);
%     xVec(2) = 100;
%     xVec(3) = 125;
    for i = 1:length(areaList)
        if ~isnan(areaList(i))
            xVec(i) = areaList(i);
        end
    end
    SA0 = pi;
    a0 = aCoeff(1);
    a1 = aCoeff(2);
    if SACur < (transSA-dx/2)*SA0%T <= .025 % resting tension -> linear increase in poly force with contact area
        timePolyForce = a0 + a1 * contactArea/xVec(2);
    elseif SACur < (transSA+dx/2)*SA0
        if isempty(coeffCell{1})
            xVec(1) = contactArea;
            areaList(1) = contactArea;
            cubicMat(1,:) = [3*xVec(1)^2, 2*xVec(1), 1, 0];
            cubicMat(2,:) = [3*xVec(3)^2, 2*xVec(3), 1, 0];
            cubicMat(3,:) = [xVec(1)^3, xVec(1)^2, xVec(1), 1];
            cubicMat(4,:) = [xVec(3)^3, xVec(3)^2, xVec(3), 1];
            cubicVec = [a1/xVec(2); (n*(1-a0-a1)/120) * ((xVec(3)-xVec(2))/120)^(n-1)...
                ; a0 + a1*xVec(1)/xVec(2); (a0+a1) + (1-a0-a1)*((xVec(3)-xVec(2))/120)^n];
            coeffCell{1} = cubicMat\cubicVec;
        end
        timePolyForce = polyval(coeffCell{1},contactArea);
        if SACur > transSA*SA0 && isnan(areaList(2))
            areaList(2) = contactArea;
        end
    else
        if isempty(coeffCell{2})
            if n == 1
                quadMat = zeros(3,3);
                quadMat(1,:) = [(xVec(2)+100)^2, xVec(2)+100, 1];
                quadMat(2,:) = [2*(xVec(2)+100), 1, 0];
                quadMat(3,:) = [(xVec(2)+140)^2, xVec(2)+140, 1];
                quadVec = [a0+a1 + (1-a0-a1)*(5/6); (1-a0-a1)/120; 1];
                coeffCell{2} = quadMat \ quadVec;
            else
                cubicMat = zeros(4,4);
                cubicMat(1,:) = [(xVec(2)+100)^3, (xVec(2)+100)^2, xVec(2)+100, 1];
                cubicMat(2,:) = [3*(xVec(2)+100)^2, 2*(xVec(2)+100), 1, 0];
                cubicMat(3,:) = [(xVec(2)+140)^3, (xVec(2)+140)^2, xVec(2)+140, 1];
                cubicMat(4,:) = [3*(xVec(2)+140)^2, 2*(xVec(2)+140), 1, 0];
                cubicVec = [a0+a1 + (1-a0-a1)*(5/6)^n; n*((5/6)^(n-1))*(1-a0-a1)/120; 1; 0];
                coeffCell{2} = cubicMat \ cubicVec;
            end
        end
        if dx == 0
            if contactArea < xVec(2)+120
                timePolyForce =  (a0+a1) + (1-a0-a1) * ((contactArea-xVec(2))/120)^n;
            else
                timePolyForce = 1;
            end
        else
            if contactArea < xVec(2)+100
                timePolyForce =  (a0+a1) + (1-a0-a1) * ((contactArea-xVec(2))/120)^n;
            elseif contactArea < xVec(2)+140
                timePolyForce = polyval(coeffCell{2},contactArea);
            else
                timePolyForce = 1;
            end
        end
    end
end

function dtVal = dtCalc(meshStack,IEN,nodesPerEl)
%     numCols = sum(zMesh <= eps);
%     numRows = length(rMesh)/numCols;
%     numRows = sum(rMesh<=eps);
%     numCols = length(rMesh)/numRows;
%     numEl = (numCols-1) * (numRows-1);
    rMesh = meshStack(:,1);
    zMesh = meshStack(:,2);
    v_rMesh = meshStack(:,3);
    v_zMesh = meshStack(:,4);
    numEl = length(IEN)/nodesPerEl;
    v_rInterp = scatteredInterpolant(rMesh,zMesh,v_rMesh);
    v_zInterp = scatteredInterpolant(rMesh,zMesh,v_zMesh);
    dtRefVals = zeros(1,numEl);
    totArea = 0;
    totDeriv = 0;
    for i = 1:numEl
%         curRow = floor((i-1) / (numCols - 1)) + 1;
%         curCol = mod(i-1, (numCols - 1)) + 1;
%         startingNode = numCols * (curRow - 1) + curCol;
%         endingNode = numCols * curRow + curCol;
%         globalIdx = [startingNode,startingNode+1,endingNode+1,endingNode];
        startIdx = nodesPerEl*(i-1)+1;
        globalIdx = IEN(startIdx:startIdx+nodesPerEl-1);
        rCur = rMesh(globalIdx);
        zCur = zMesh(globalIdx);
        gaussPts = [-1/sqrt(3), 1/sqrt(3)];
        rC = 0;
        zC = 0;
        area = 0;
        testDeriv = 0;
        for n = 1:length(gaussPts)
            for m = 1:length(gaussPts)
                curEps = gaussPts(n);
                curEta = gaussPts(m);
                N_eps = [ -(1 + curEta) / 4, (1 + curEta) / 4, (1 - curEta) / 4, -(1 - curEta) / 4 ]';
                N_eta = [ (1 - curEps) / 4, (1 + curEps) / 4, -(1 + curEps) / 4, -(1 - curEps) / 4 ]';
                N_cur = [ 0.25 * (1 - curEps) * (1 + curEta), 0.25 * (1 + curEps) * (1 + curEta),...
                    0.25 * (1 + curEps) * (1 - curEta), 0.25 * (1 - curEps) * (1 - curEta) ]';
                if nodesPerEl == 3
                    N_eps = [N_eps(1:2); N_eps(3)+N_eps(4)];
                    N_eta = [N_eta(1:2); N_eta(3)+N_eta(4)];
                    N_cur = [N_cur(1:2); N_cur(3)+N_cur(4)];
                end
                x_eps = sum(N_eps .* rCur);
                x_eta = sum(N_eta .* rCur);
                y_eps = sum(N_eps .* zCur);
                y_eta = sum(N_eta .* zCur);
                jVal = x_eps * y_eta - x_eta * y_eps;
                jMat = (1/jVal) * [y_eta, -x_eta; -y_eps, x_eps];
                NMat = horzcat(N_eps, N_eta);
                NMat = NMat * jMat; %  4x2, each row a is [dN_a/dx, dN_a/dy]
                testDeriv = testDeriv + jVal*(sum(NMat(:,1).*sin(pi*rCur/.5)) + sum(NMat(:,2).*sin(pi*zCur)));
                rC = rC + jVal*sum(N_cur .* rCur); 
                zC = zC + jVal*sum(N_cur .* zCur);
                area = area + jVal;
            end
        end
        rC = rC/area;
        zC = zC/area;
%         if 1
%             curPolygon = polyshape(rCur,zCur);
%             [rCRef,zCRef] = centroid(curPolygon);
%             areaRef = polyarea(rCur,zCur);
%         end
        v_rCur = v_rInterp(rC,zC);
        v_zCur = v_zInterp(rC,zC);
        rDim = max(rCur) - min(rCur);
        zDim = max(zCur) - min(zCur);
        dtRefVals(i) = min([abs(rDim/v_rCur), abs(zDim/v_zCur)]);
%         vMag = sqrt(v_rCur^2 + v_zCur^2);
%         charLength = sqrt(area/2);
%         dtRefVals(i) = charLength/vMag;
        totArea = totArea + area;
        totDeriv = totDeriv + testDeriv;
    end
    dtVal = min(dtRefVals);
%     assessLogic = ~isnan(dtRefVals);
%     dtVal = median(dtRefVals(assessLogic));
end

function [r,z,prevSurfLogic,errorFlag] = discreteThresh(r,z,NLoc,prevSurfLogic,T,SA,receptorDensity,diffuseRange)
    notStuck = ~prevSurfLogic & z <= eps;
    startSurf = find(prevSurfLogic,1,'first');
    foldBack = notStuck & r < r(startSurf);
    r = r(~foldBack);
    z = z(~foldBack);
    prevSurfLogic = prevSurfLogic(~foldBack);
    KbT = 4.11e-18 / 1e-12; %KbT for T = 25 deg C, in units of nN * microns
    tau = T; % tension in nN/micron (same as mN/m)
    bendMod = 1e-3; % bending modulus in nN * microns
    threshDist = sqrt((1/(4*pi*tau/KbT)) * ...
        log(1 + ((tau*SA*(8.5^2))/(pi^2*bendMod))));
    threshDist = threshDist/8.5;
    startAttachIdx = find(z <= threshDist & ~prevSurfLogic,1,'first');
    if isempty(startAttachIdx)
        testSpan = startSurf-1:startSurf;
    else
        testSpan = startAttachIdx-1:startSurf; % length of membrane close to surface
    end
    rSpan = r(testSpan);
    zSpan = z(testSpan);
    zSpan(zSpan<0) = 1e-12;
    if ~isinf(receptorDensity)
        zSpan(zSpan<0) = 1e-12;
        sParam = zeros(1,length(rSpan));
        for i = 2:length(sParam)
            sParam(i) = sParam(i-1) + sqrt((rSpan(i)-rSpan(i-1))^2 + (zSpan(i)-zSpan(i-1))^2);
        end
        receptorLoc = sParam(end):-(1/sqrt(receptorDensity)):sParam(1);
        if length(receptorLoc) == 1
            receptorLoc = [];
        else
            receptorLoc = receptorLoc(2:end);
            rReceptor = interp1(sParam,rSpan,receptorLoc);
            zReceptor = interp1(sParam,zSpan,receptorLoc);
            zReceptor(zReceptor < 0) = 0;
        end
    end
    errorFlag = false;
    if ~isempty(rSpan)
        NTestLogic = NLoc > rSpan(end)+eps & NLoc <= rSpan(1);
        NTest = NLoc(NTestLogic); % possible new binding sites
%         bindLogic = bindLogic(NTestLogic);
        rSpanNew = rSpan;
        zSpanNew = zSpan;
        contactLogic = false;

        for i = 1:length(NTest)
            curDist = inf;
            if isinf(receptorDensity)
                for j = 2:length(rSpan) % for each ligand, test distance to cell membrane (segment by segment)
                    rSeg = rSpan(j-1:j);
                    zSeg = zSpan(j-1:j);
                    r1 = rSeg(1);
                    z1 = zSeg(1);
                    dr = rSeg(2) - rSeg(1);
                    dz = zSeg(2) - zSeg(1);
                    % find minimum distance between current ligand and
                    % this membrane segment
                    sVal = ((NTest(i)-r1)*dr - z1*dz) / (dz^2 + dr^2); % rel location of min
                    % (s is a parameter where s=0 is start of current
                    % segment and s=1 is the end of the segment)
                    if sVal > 1 % then min is beyond segment
                        sVal = 1;
                    elseif sVal < 0 % then min is before segment
                        sVal = 0;
                    end
                    testDist = sqrt((r1+sVal*dr-NTest(i))^2 + (z1+sVal*dz)^2);
                    if testDist < curDist % then this is the new minimum distance
                        curDist = testDist;
                        [~,distIdx] = min(abs(rSpanNew-rSeg(1)));
                    end
%                     if any(zSeg <= 0 & rSeg > NTest(i))
%                         contactLogic = true; % membrane comes in contact beyond current binding site -> consider bound at site
%                     end
%                     if zSeg(1) > 0 && zSeg(2) < 0
%                         try rZero = interp1(zSeg,rSeg,0);
%                         catch
%                             errorFlag = true;
%                             break
%                         end
%                         if rZero > NTest(i)
%                             contactLogic = true; % membrane comes in contact beyond current binding site -> consider bound at site
%                         end
%                     end
                end
            else
                for j = 1:length(receptorLoc)
%                     if diffuseRange ~= 0
                        curRange = [ max([receptorLoc(j)-diffuseRange,sParam(1)]), ...
                            min([receptorLoc(j)+diffuseRange, sParam(end)]) ];
                        rSeg = interp1(sParam,rSpan,curRange);
                        zSeg = interp1(sParam,zSpan,curRange);
                        if zReceptor == 0
                            rSeg = [max([0, rReceptor(j)-diffuseRange]), rReceptor(j)+diffuseRange];
                            zSeg = [0,0];
                        end
                        r1 = rSeg(1);
                        z1 = zSeg(1);
                        dr = rSeg(2) - rSeg(1);
                        dz = zSeg(2) - zSeg(1);
                        % find minimum distance between current ligand and
                        % this membrane segment
                        sVal = ((NTest(i)-r1)*dr - z1*dz) / (dz^2 + dr^2); % rel location of min
                        % (s is a parameter where s=0 is start of current
                        % segment and s=1 is the end of the segment)
                        if sVal > 1 % then min is beyond segment
                            sVal = 1;
                        elseif sVal < 0 % then min is before segment
                            sVal = 0;
                        end
                        testDist = sqrt((r1+sVal*dr-NTest(i))^2 + (z1+sVal*dz)^2);
                        if testDist < curDist
                            curDist = testDist;
                            distIdx = find(sParam <= curRange(1),1,'last');
                        end
%                     else
%                         testDist = sqrt((rReceptor(j)-NTest(i))^2 + (zReceptor(j))^2);
%                         if testDist < curDist
%                             curDist = testDist;
%                             distIdx = find(sParam <= receptorLoc(j),1,'last');
%                         end
%                     end
                end
            end
            if (curDist < threshDist) % new binding occurs! yay!
                rSpanNew = [rSpanNew(1:distIdx),NTest(i),rSpanNew(distIdx+1:end)];
                zSpanNew = [zSpanNew(1:distIdx),0,zSpanNew(distIdx+1:end)];
            end
        end
        % add on bound points to vectors
        r = [r(1:testSpan(1)-1),rSpanNew,r(startSurf+1:end)];
        z = [z(1:testSpan(1)-1),zSpanNew,z(startSurf+1:end)];
        surfLogicNewSpan = zSpanNew == 0;
        prevSurfLogic = [prevSurfLogic(1:testSpan(1)-1),surfLogicNewSpan,prevSurfLogic(startSurf+1:end)];
        notStuck = ~prevSurfLogic & z <= eps; % "on surface" but unbound -> does not stick
        z(notStuck) = 1e-12;
    end
end