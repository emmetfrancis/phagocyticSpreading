i = 50;
[rTest,zTest,elDissip,totDissip,splitDissip] = dissipationCalc(rStored{i},zStored{i},v_rStored{i},v_zStored{i});
surfLogic = zStored{i}<=eps;
numCols = sum(surfLogic) - 1;
numRows = length(rTest)/numCols;
rContact = rTest(end);
rTest = [rTest,(rContact/(numCols-1))*(0:numCols-1)];
zTest = [zTest,zeros(1,numCols)];
elDissip = [elDissip,zeros(1,numCols)];
numRows = numRows+1;
rMesh = reshape(rTest,[numCols,numRows]);
zMesh = reshape(zTest,[numCols,numRows]);
rMesh = rMesh';
zMesh = zMesh';
elDissipMat = reshape(abs(elDissip),[numCols,numRows]);
elDissipMat = elDissipMat';
figure
surf(rMesh,zMesh,elDissipMat,elDissipMat,'FaceColor','interp','EdgeColor','none')
caxis([0 7.2464e5])
xlim([.15 .32])
ylim([0 .08])
daspect([1 1 1])
view([0 90])
prettyGraph
grid off
hold on
[rBound,zBound] = meshToBound(rTest,zTest);
plot(rBound,zBound,'k','LineWidth',2)

%%
figure
for i = 1:6
    subplot(2,3,i)
    curSplit = splitDissip{i};
%     if i==4
%         curSplit = splitDissip{i} + splitDissip{i+1} + splitDissip{i+2};
%     end
%     curInterp = scatteredInterpolant(rTest',zTest',curSplit');
%     dissipCur = curInterp(rCur,zCur);
%     FEMPlot_p(rCur,zCur,dissipCur,numRows)
%     title(sprintf('Dissip portion %d',i))
%     caxis([0 max(elDissip)/10])
    curSplit = [curSplit,zeros(1,numCols)];
    curSplitMat = reshape(abs(curSplit),[numCols,numRows]);
    curSplitMat = curSplitMat';
    surf(rMesh,zMesh,curSplitMat,curSplitMat,'FaceColor','interp','EdgeColor','none')
    caxis([0 max(elDissip)/10])
    zlim([0 max(elDissip)])
    xlim([.15 .32])
    ylim([0 .08])
    daspect([1 1 1])
    view([0 90])
    prettyGraph
    title(sprintf('case %d',i))
    grid off
end

%% plot z dissipation profile
curLogic = rTest > rc(i) & rTest < rc(i)+.05;
[~,maxIdx] = max(elDissip);
zInterp = 0.01:.005:0.2;
rInterp = rTest(maxIdx)*ones(size(zInterp));
dissipInterp = scatteredInterpolant(rTest',zTest',elDissip');
zProfileDissip = dissipInterp(rInterp,zInterp);
figure
% loglog(zTest(curLogic),elDissip(curLogic)./(rc(i)*radGrowth(i)),'*')
loglog(zInterp,zProfileDissip,'*')
hold on
% curLogic = curLogic & zTest > .01 & zTest < .05;
% figure
% plot(rTest(curLogic)-rc(i),sqrt(elDissip(curLogic))./(rc(i)*radGrowth(i)),'*')
% hold on

%% calc dissip
dissipOverTime = zeros(length(v_rStored),1);
% rWidth = zeros(length(v_rStored),1);
zMaxVec = zeros(length(v_rStored),1);
zShearMax = zeros(length(v_rStored),1);
mixedDissipTot = zeros(length(v_rStored),1);
for i = 1:length(v_rStored)
    [rTest,zTest,elDissipCur,dissipOverTime(i),splitDissip] = dissipationCalc(rStored{i},zStored{i},v_rStored{i},v_zStored{i});
%     sampleLogic = rTest > rc(i) & rTest < rc(i)+0.02 & zTest < 0.02;
%     zShearMax(i) = mean(splitDissip{6}(sampleLogic));
%     mixedDissipTot(i) = sum(splitDissip{5}(sampleLogic));
%     curLogic = rTest > rc(i)-.1 & rTest < rc(i)+.1;
%     [rMax,~] = max(elDissipCur(curLogic));
%     [~,rMaxIdx] = min(abs(elDissipCur-rMax));
%     rCenter = rTest(rMaxIdx);
%     rOuter = max(rTest(elDissipCur > rMax*.1));
%     rInner = min(rTest(elDissipCur > rMax*.1));
%     rCenterVec(i) = rCenter;
%     rWidth(i) = rOuter - rInner;
%     curLogic = curLogic & zTest > .01;
%     
%     zMaxVec(i) = max(elDissipCur(curLogic));
% %     [~,zMaxIdx] = min(abs(elDissipCur-zMax));
% %     zMaxVec(i) = zTest(zMaxIdx);
%     zDissipTot(i) = sum(elDissipCur(curLogic));
end

% timeVec = timeVec/170;
SA = zeros(length(rStored),1);
vol = zeros(length(rStored),1);
rc = zeros(length(rStored),1);
tension = zeros(length(rStored),1);
pEnergy = zeros(length(rStored),1);
SA0 = pi;
dx = 0.26; % determines smoothing of tension curve
transSA = 1.26;
for i = 1:length(rStored)-1
    [rBound,zBound] = meshToBound(rStored{i},zStored{i});
    [SA(i),vol(i)] = SAVolCalc(rBound,zBound,'trapz');
    surfLogic = zBound <= eps;
    rc(i) = max(rBound(surfLogic));
    if SA(i) < (transSA-dx/2)*SA0
        tension(i) = .01 + .16*(SA(i)-SA0)/SA0;
    elseif SA(i) < (transSA+dx/2)*SA0
        SAArg = (SA(i)-SA0)/SA0;
        tension(i) =  .7405*SAArg^2 - 0.0325*SAArg + 0.0225;
    else
        tension(i) = .0516 + .54506*(SA(i)-transSA*SA0)/SA0;
    end
    pEnergy(i) = pressureEnergy(rStored{i},zStored{i},pStored{i});
end
tension = tension/.01;
figure
loglog(timeVec(1:length(v_rStored)),dissipOverTime)
hold on
tensionEnergy = smooth(diff(SA.*tension)./diff(timeVec'),31);
pEnergyIncrease = smooth(diff(pEnergy)./diff(timeVec'),31);
loglog(timeVec(1:end-1),tensionEnergy)
loglog(timeVec(1:end-1),dissipOverTime(1:length(tensionEnergy))+tensionEnergy)

%%
energyGain = smooth(diff(pi*rc.^2)./diff(timeVec'),31);
radGrowth = smooth(diff(rc)./diff(timeVec'),31);
allEnergy = dissipOverTime(1:length(tensionEnergy))+tensionEnergy+pEnergyIncrease;
figure
yyaxis left
loglog(timeVec(1:length(tensionEnergy)),allEnergy)
yyaxis right
loglog(timeVec(1:end-1),energyGain)

%% relate contact area growth to velocities
figure
dz = .01;
vrFront = zeros(length(v_rStored),1);
vzFront = zeros(length(v_zStored),1);
for i = 1:length(v_rStored)
    curLogic = zStored{i} > dz & zStored{i} < 2*dz & rStored{i} > rc(i)-dz & rStored{i} < rc(i)+dz;
    [~,sampleIdx] = min(abs(rStored{i}(curLogic)-rc(i)));
    vrFront(i) = mean(v_rStored{i}(sampleIdx));
    vzFront(i) = mean(v_zStored{i}(sampleIdx));
end
% semilogy(timeVec(1:length(v_rStored)),vrFront)
% hold on
% semilogy(timeVec(1:length(v_rStored)),abs(vzFront))
% semilogy(timeVec(1:end-1),radGrowth)
loglog(rc(1:length(v_zStored)),abs(vzFront)./radGrowth,'*')

%% plot v_z profile
figure
dz = .01;
for i = 1:10:100%length(v_zStored)
    rCur = rStored{i};
    zCur = zStored{i};
    v_zCur = v_zStored{i};
%     rInterp = (rc(i)-.1:.001:rc(i)+.2)';
%     zInterp = .01*ones(size(rInterp));
    zInterp = (.01:.005:.5);
    rInterp = 1.05*rc(i)*ones(size(zInterp));
    curInterp = scatteredInterpolant(rCur,zCur,v_zCur);
    v_zInterp = curInterp(rInterp,zInterp);
    loglog(zInterp,-v_zInterp,'*')
%     loglog(rInterp-rc(i),-v_zInterp,'*')
%     v_zGrad = smooth(diff(v_zInterp)./(diff(rInterp)),5);
%     loglog((rInterp(1:end-1)-rc(i))./rc(i),-v_zGrad/abs(mean(v_zInterp(101:111))),'*')
    hold on
end

%%
rc = 1;
Fz = 1;
mu = 1;
r = .99*rc;
z = 0:.01:.5;
[r,z] = meshgrid(r,z);
R = sqrt((r-rc).^2 + z.^2);
vr = ((Fz.*(rc-r))./(4*mu)).*2.*z.*R.^(-3);
vz = -(1./r)*(Fz/(4*mu)) .* ((r-rc)./R + r.*((r-rc).^2+z.^2-(r-rc)).*R.^(-3));
% figure
% quiver(r,z,vr,vz)
figure
plot(z,vr,'*')

%% analytical solution
rcVals = 0.11:.1:5;
dissipVals = zeros(size(rcVals));
for i = 1:length(rcVals)
    rc = rcVals(i);
    syms r z
    R = sqrt((r-rc).^2 + z.^2);
    psi = rc.*(r-rc).^2 ./ R;
    vz = (1/r) * diff(psi,r);
    dissip = (diff(vz,r).^2);
    rVec = rc-.1:.01:rc+.1;
    zVec = 0.01:.001:.03;
    [rMesh,zMesh] = meshgrid(rVec,zVec);
    rMesh = rMesh(:);
    zMesh = zMesh(:);
    dissipMesh = zeros(size(rMesh));
    for j = 1:length(dissipMesh)
        dissipMesh(j) = subs(dissip,[r,z],[rMesh(j),zMesh(j)]);
    end
    rMesh = reshape(rMesh,[length(zVec),length(rVec)]);
    zMesh = reshape(zMesh,[length(zVec),length(rVec)]);
    dissipMesh = reshape(dissipMesh,[length(zVec),length(rVec)]);
%     dissipIntArg = subs(dissip.*r,r,rc);
%     dissipInt1 = int(dissip.*r,r,[rc-.01 rc+.01]);
%     dissipInt2 = int(dissipInt2,z,[.01,10]);
%     dissipVals(i) = double(dissipInt2);
    dissipVals(i) = trapz(zVec,trapz(rVec,rMesh.*dissipMesh,2));
end

%% another analytical test
syms phi_c r
% phi_c = 0:pi/100:2*pi;
rVec = .8:.1:1.2;
zVec = .1:.2:10;
% phiVec = 0:pi/5:2*pi;
phiVec = 0;
r_c = .5;
% [r,z,phi] = meshgrid(rVec,zVec,phiVec);
% r = r(:);
z = zVec(:);
phi = phiVec(:);
% z = z(:);
% phi = phi(:);
vzAll = zeros(size(phi));
% vzAllPsi = zeros(size(phi));
dissipVec = zeros(size(z));
for i = 1:length(z)
%     rCur = r(i);
    rCur = r;
    zCur = z(i);
    phiCur = phi;%(i);
    rTrans = sqrt((rCur.*cos(phiCur)-r_c.*cos(phi_c)).^2 + (rCur.*sin(phiCur)-r_c.*sin(phi_c)).^2);
    R = sqrt(rTrans.^2 + zCur.^2);
    % R = sqrt(((r-r_c)*cos(phi_c)).^2 + z.^2);
    vz_Stokeslet = (1./R + zCur.^2./R.^3);
%     psi = rTrans.^2./R;
%     psiInt = int(psi,phi_c,[0 2*pi]);
%     vzFromPsi = (1./r)*diff(psiInt,r);
    vzInt = int(vz_Stokeslet,phi_c,[0 2*pi]);
    vzAll(i) = double(subs(vzInt,r,r_c));
%     vzAllPsi(i) = double(subs(vzFromPsi,r,0.9));
%     vzAll(i) = trapz(phi_c,vz_Stokeslet);
    dissip = diff(vzInt,r);
    dissipVec(i) = double(subs(dissip,r,r_c));
end

%% quick geometry calc
rc = 0.01:.001:.4;
r0 = .05/8.5;
% R = 0.5;
theta_c = asin(rc./R);
R = (4*0.5^3./((1+cos(theta_c)).^2 .* (2-cos(theta_c)))).^(1/3);
threshAng = acos(cos(theta_c) - r0./R);
dr = R.*(sin(threshAng)-sin(theta_c));
figure
loglog(rc,dr)

%% prettier log log plot
xVals = rc(1:end-1).^2 .* radGrowth;
xVals = xVals*(8.5e-6^2) * (1e-5/200) * 1e18; % convert to um^3 / s
xEdges = logspace(-2,0.7,50);
xSample = (xEdges(1:end-1) + xEdges(2:end))/2;
dissipUnits = dissipOverTime * 1e18 * (200*((1e-5/200)/8.5e-6)^2 * (8.5e-6)^3);
dissipSample = cell(size(xSample));
meanDissipSample = zeros(size(xSample));
stdDissipSample = zeros(size(xSample));
for i = 1:length(xSample)
    curLogic = find(xVals > xEdges(i) & xVals < xEdges(i+1));
    dissipSample{i} = dissipUnits(curLogic);
    meanDissipSample(i) = mean(dissipSample{i});
    stdDissipSample(i) = std(dissipSample{i});
end
% fitLogic = xSample > 0.5 & xSample < 3;
% logFit = polyfit(log(xSample(fitLogic)),log(meanDissipSample(fitLogic)),1);
figure
% loglog(xSample,exp(polyval(logFit,log(xSample))))
loglog(xSample,650*xSample.^2)
xlim([.1 3])
ylim([10 1e4])
hold on
errorbar(xSample,meanDissipSample,stdDissipSample,'LineStyle','none')

% scatter(xVals,dissipOverTime,[],timeVec(1:end-1))
% colorbar

%% just doing things
for i = 1:length(v_rStored)
    [rTest,zTest,elDissip,totDissip,splitDissip] = dissipationCalc(rStored{i},zStored{i},v_rStored{i},v_zStored{i});
    zShearTot(i) = sum(splitDissip{6});
    allShearTot(i) = sum(splitDissip{4}+splitDissip{5}+splitDissip{6});
end
    