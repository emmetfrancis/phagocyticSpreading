function [vx,vy,KStruct,fFull,fBCFull] = FEMStokesVelocityCalc(meshStruct,pVec,displayLogic,KStruct,fBCFull)

P = meshStruct.P;
IEN = meshStruct.IEN;
x = meshStruct.x;
y = meshStruct.y;
gVec = meshStruct.gVec;
hVec = meshStruct.hVec;
gNum = meshStruct.gNum;
axisymm = meshStruct.axisymm;
if isfield(meshStruct,'numRows') % then quad elements
    numRows = meshStruct.numRows;
    numCols = meshStruct.numCols;
    nodesPerEl = 4;
    numEl = (numRows-1)*(numCols-1);
    dof = 2*numRows * numCols - gNum;
else % then triangles
    nodesPerEl = 3;
    numEl = meshStruct.numEl;
    dof = 2*max(IEN) - gNum;
end

if isempty(KStruct)
    K = zeros(dof, dof);
    fBCFull = zeros(dof,1);
    KCalc = false;
else
    K = KStruct{1};
    Kinv = KStruct{2};
    KCalc = true;
end
xVec = zeros(max(IEN),1);
yVec = zeros(max(IEN),1);
for i = 1:length(IEN)
    xVec(IEN(i)) = x(i);
    yVec(IEN(i)) = y(i);
end
% figure
u0 = zeros(2*length(xVec),1);
% FEMPlot(xVec,yVec,u0,numRows)
% hold on
T0 = 1;
% T0 = T0*1000;
f = zeros(dof,1);
if isempty(pVec) % then initialize with uniform pressure
%     [coeffFile,path]= uigetfile('*.mat','Choose file for coefficient values');
%     coeffLoad = load(fullfile(path,coeffFile));
%     coeffVals = coeffLoad.coeffStored{2};
%     [~,~,pVec] = meshFlowCalc(xVec,yVec,coeffVals,'odd');
    pVec = 4.015*ones(size(xVec));
%     pVec(abs(y)>L/2) = 2/Rp;
%     pVec(abs(y)<=L/2) = 1/Rp;
end
% pVec = (T0/8.5)*pVec;
p = zeros(length(IEN),1);
for i = 1:length(IEN)
    p(i) = pVec(IEN(i));
end
% p = zeros(length(p),1);
% p = p - min(p);
% p = (8.5/T0)*ones(size(p));
% p(round(length(p)/2)+1:end) = 0;
% p = 0*p;
mu = 1;
% hVec = hVec*0; %tensionless case
% p = abs(y);
gVecFull = zeros(length(f),1);
if length(hVec) == length(x)
    tensionOnly = true;
elseif length(hVec) == length(P)
    tensionOnly = false;
else
    error('length hVec does not match x or P!')
end
for i = 1:numEl
    PIdx = (i-1) * 2*nodesPerEl + 1;
    xIdx = (i-1) * nodesPerEl + 1;
    xCur = x(xIdx:xIdx+(nodesPerEl-1));
    yCur = y(xIdx:xIdx+(nodesPerEl-1));
    if ~KCalc
        K_el = elStiffness(mu, xCur, yCur, axisymm);
        K = KAssemble(P(PIdx:PIdx+(2*nodesPerEl-1)), K_el, K);
        if tensionOnly
            fBC = elForceBC(hVec(xIdx:xIdx+(nodesPerEl-1)), T0, xCur, yCur, axisymm);
        else
            fBC = elForceBC(hVec(PIdx:2:PIdx+(2*nodesPerEl-2)),...
                hVec(PIdx+1:2:PIdx+(2*nodesPerEl-1)),xCur,yCur,axisymm);
        end
        fBCFull = fAssemble(P(PIdx:PIdx+(2*nodesPerEl-1)), fBC, fBCFull);
        g_fVec = K_el * gVec(PIdx:PIdx+(2*nodesPerEl-1));
    else
        g_fVec = zeros(2*nodesPerEl,1);
    end
%     gVecFull = fAssemble(P(PIdx:PIdx+7),gVec(PIdx:PIdx+7),gVecFull);
    fEl = elForce(p(xIdx:xIdx+(nodesPerEl-1)), xCur, yCur, axisymm, g_fVec);
%     fEl = fEl + fBC;
    f = fAssemble(P(PIdx:PIdx+(2*nodesPerEl-1)), fEl, f);
end
% fBCFull(1:2:end) = abs(fBCFull(1:2:end));
% halfLength = round(length(fBCFull)/2);
% fBCFull(halfLength+2:2:end) = -fBCFull(halfLength+2:2:end);
% fFull = -2*f + fBCFull;
fFull = -f + fBCFull;
% fFull = fFull/4;
% tic
% uSolve = pinv(K) * fFull;
% uSolve = K \ fFull;
if ~KCalc
%     K = sparse(K);
    Kinv = inv(K);
    if any(isnan(Kinv(:))) || any(isinf(Kinv(:)))
        Kinv = pinv(K);
    end
% else
%     Kinv = inv(K);
end
uSolve = Kinv * fFull;
% toc
% pShift = 0;
% pCorr = 1;
% meanResid = 1;
% while abs(meanResid) > .1
%     fFull = -2*pCorr*f + fBCFull;
%     fFull = fFull*(.025/.2);
%     uVec = pinv(K) * fFull;
%     vx = uVec(1:2:end);
%     vy = uVec(2:2:end);
%     dvResid = checkIncompressible(meshStruct,vx,vy);
%     meanResid = mean(dvResid);
%     pShift = 
% assemble global u vec
uVec = u0;
% yLogic = false(length(uVec),1);
for i = 1:length(uSolve)
    globalIdx = find(P == i,1,'first');
    nodeNum = IEN(ceil(globalIdx/2));
    xyNum = mod(globalIdx,2);
    if xyNum == 1
        uVec(2*nodeNum-xyNum) = uSolve(i);
    else
        uVec(2*nodeNum) = uSolve(i);
%         yLogic(2*nodeNum) = true;
    end
end
% correct for EBC
% uVec(yLogic) = uVec(yLogic) - mean(uVec(yLogic));
if displayLogic
    FEMPlot(xVec,yVec,uVec,numRows)
end
vx = uVec(1:2:end);
vy = uVec(2:2:end);
KStruct{1} = K;
KStruct{2} = Kinv;
end


function K_el = elStiffness(mu, x, y, axisymm) 
    numNodes = length(x);
    K_el = zeros(numNodes*2, numNodes*2);
%     gaussPts = [-sqrt(3/5),0,sqrt(3/5)];
%     gaussWts = [5/9,8/9,5/9];
    gaussPts = [-1/sqrt(3), 1/sqrt(3)];
    gaussWts = [1,1];
    for n = 1:length(gaussPts)
        for m = 1:length(gaussPts)
            curEps = gaussPts(n);
            curEta = gaussPts(m);
            N_eps = [ -(1 + curEta) / 4, (1 + curEta) / 4, (1 - curEta) / 4, -(1 - curEta) / 4 ]';
            N_eta = [ (1 - curEps) / 4, (1 + curEps) / 4, -(1 + curEps) / 4, -(1 - curEps) / 4 ]';
            N_cur = [ 0.25 * (1 - curEps) * (1 + curEta), 0.25 * (1 + curEps) * (1 + curEta),...
                0.25 * (1 + curEps) * (1 - curEta), 0.25 * (1 - curEps) * (1 - curEta) ]';
            if numNodes == 3
                N_eps = [N_eps(1:2); N_eps(3)+N_eps(4)];
                N_eta = [N_eta(1:2); N_eta(3)+N_eta(4)];
                N_cur = [N_cur(1:2); N_cur(3)+N_cur(4)];
            end
            x_eps = sum(N_eps .* x);
            x_eta = sum(N_eta .* x);
            y_eps = sum(N_eps .* y);
            y_eta = sum(N_eta .* y);
            jVal = x_eps * y_eta - x_eta * y_eps;
            jMat = (1/jVal) * [y_eta, -x_eta; -y_eps, x_eps];
            NMat = horzcat(N_eps, N_eta);
            NMat = NMat * jMat; %  4x2, each row a is [dN_a/dx, dN_a/dy]
            gWt = gaussWts(n)*gaussWts(m);
            for i = 1:numNodes
                for j = i:numNodes
                    idx1 = 2*i-1:2*i;
                    idx2 = 2*j-1:2*j;
                    elSubMat = zeros(2,2);
                    elSubMat(1,2) = (mu/4)*NMat(i,2)*NMat(j,1);
                    elSubMat(2,1) = (mu/4)*NMat(i,1)*NMat(j,2);
                    elSubMat(2,2) = (mu/2)*NMat(i,2)*NMat(j,2) + (mu/4)*NMat(i,1)*NMat(j,1);
                    if axisymm
                        xCur = sum(N_cur .* x);
                        elSubMat(1,1) = (mu/2)*NMat(i,1)*NMat(j,1) + (mu/4)*NMat(i,2)*NMat(j,2) +...
                            (mu/2)*N_cur(i)*N_cur(j)/xCur^2;
                        elSubMat = 4*elSubMat;
                        K_el(idx1,idx2) = K_el(idx1,idx2) + gWt*jVal * xCur * elSubMat;
                        if i ~= j
                            K_el(idx2,idx1) = K_el(idx2,idx1) + gWt*jVal * xCur * elSubMat';
                        end
                    else
                        elSubMat(1,1) = (mu/2)*NMat(i,1)*NMat(j,1) + (mu/4)*NMat(i,2)*NMat(j,2);
                        K_el(idx1,idx2) = K_el(idx1,idx2) + gWt*jVal * elSubMat;
                        if i ~= j
                            K_el(idx2,idx1) = K_el(idx2,idx1) + gWt*jVal * elSubMat';
                        end
                    end
                end
            end
        end
    end
end

function fVec = elForce(p, x, y, axisymm, g_fVec)
    numNodes = length(x);
    fVec = zeros(2*numNodes,1);
%     gaussPts = [-sqrt(3/5),0,sqrt(3/5)];
%     gaussWts = [5/9,8/9,5/9];
    gaussPts = [-1/sqrt(3), 1/sqrt(3)];
    gaussWts = [1,1];
    for n = 1:length(gaussPts)
        for m = 1:length(gaussPts)
            curEps = gaussPts(n);
            curEta = gaussPts(m);
            N_eps = [ -(1 + curEta) / 4, (1 + curEta) / 4, (1 - curEta) / 4, -(1 - curEta) / 4 ]';
            N_eta = [ (1 - curEps) / 4, (1 + curEps) / 4, -(1 + curEps) / 4, -(1 - curEps) / 4 ]';
            N_cur = [ 0.25 * (1 - curEps) * (1 + curEta), 0.25 * (1 + curEps) * (1 + curEta),...
                0.25 * (1 + curEps) * (1 - curEta), 0.25 * (1 - curEps) * (1 - curEta) ]';
            if numNodes == 3
                N_eps = [N_eps(1:2); N_eps(3)+N_eps(4)];
                N_eta = [N_eta(1:2); N_eta(3)+N_eta(4)];
                N_cur = [N_cur(1:2); N_cur(3)+N_cur(4)];
            end
            x_eps = sum(N_eps .* x);
            x_eta = sum(N_eta .* x);
            y_eps = sum(N_eps .* y);
            y_eta = sum(N_eta .* y);
            jVal = x_eps * y_eta - x_eta * y_eps;
            jMat = (1/jVal) * [y_eta, -x_eta; -y_eps, x_eps];
            NMat = horzcat(N_eps, N_eta);
            NMat = NMat * jMat; %  4x2, each row a is [dN_a/dx, dN_a/dy]
            % need to calculate p value from current eta and xi values
            pVal = sum(N_cur .* p);
            gWt = gaussWts(n)*gaussWts(m);
            if axisymm
                xCur = sum(N_cur .* x);
                NMat(:,1) = NMat(:,1) + N_cur/xCur;
                NMatTrans = NMat';
                fCur = pVal*NMatTrans(:);
                fVec = fVec - gWt*fCur*jVal*xCur;
            else
                NMatTrans = NMat';
                fCur = pVal*NMatTrans(:);
                fVec = fVec - gWt*fCur*jVal;
            end
        end
    end
    % subtract contribution from EBC
    fVec = fVec - g_fVec;
end


function fVec = elForceBC(xMag,yMag,x,y,axisymm)
% need to do a 1d integral along each side with an NBC
    numNodes = length(x);
    fVec = zeros(numNodes*2,1);
    epsVec = [-1, 1, 1, -1];
    etaVec = [1, 1, -1, -1];
    gaussPts = [-1/sqrt(3), 1/sqrt(3)];
    for n = 1:numNodes
        if xMag(n) ~= 0
            if n == numNodes
                nextIdx = 1;
            else
                nextIdx = n + 1;
            end
            for m = 1:2
                if epsVec(n) == epsVec(nextIdx)
                    curEps = epsVec(n);
                    curEta = gaussPts(m);
                else % then eta is constant
                    curEps = gaussPts(m);
                    curEta = etaVec(n);
                end
                dx = x(nextIdx) - x(n);
                dy = y(nextIdx) - y(n);
                dL = sqrt(dx^2 + dy^2);
                eps1D = gaussPts(m);
                N1DVal = 0.5*[1-eps1D;1+eps1D];
                j1D = dL/2;
                n_x = -dy/dL;
                n_y = dx/dL;
                if n_x < 0 
%                     error('x normal cannot be negative')
                    n_x = 0;
                end
                if length(yMag) == 1 % then this is tension only
                    h_x = -yMag * xMag(n) * n_x; % here xMag is curv
                    h_y = -yMag * xMag(n) * n_y;
                else % includes other cases
                    h_x = xMag(n);
                    h_y = yMag(n);
                end
                if axisymm
                    N_cur = [ 0.25 * (1 - curEps) * (1 + curEta), 0.25 * (1 + curEps) * (1 + curEta),...
                        0.25 * (1 + curEps) * (1 - curEta), 0.25 * (1 - curEps) * (1 - curEta) ]';
                    if numNodes == 3
                        N_cur = [N_cur(1:2); N_cur(3)+N_cur(4)];
                    end
                    xCur = sum(N_cur .* x);
%                     xCur = 1;
                    fVec(2*n-1) = fVec(2*n-1) + j1D * h_x * xCur * N1DVal(1);
                    fVec(2*n) = fVec(2*n) + j1D * h_y * xCur * N1DVal(1);
                    fVec(2*nextIdx-1) = fVec(2*nextIdx-1) + j1D * h_x * xCur * N1DVal(2);
                    fVec(2*nextIdx) = fVec(2*nextIdx) + j1D * h_y * xCur * N1DVal(2);
                else
                    fVec(2*n-1) = fVec(2*n-1) + j1D * h_x * N1DVal(1);
                    fVec(2*n) = fVec(2*n) + j1D * h_y * N1DVal(1);
                    fVec(2*nextIdx-1) = fVec(2*nextIdx-1) + j1D * h_x * N1DVal(2);
                    fVec(2*nextIdx) = fVec(2*nextIdx) + j1D * h_y * N1DVal(2);
                end
            end
        end
    end
end


function K = KAssemble(P, K_el, K)
    for n = 1:length(P)
        if P(n) > 0
            for m = 1:length(P)
                if P(m) > 0
                    K(P(n), P(m)) = K(P(n), P(m)) + K_el(n, m);
                end
            end
        end
    end
end


function f = fAssemble(P, f_el, f)
    for n = 1:length(P)
        if P(n) > 0
            f(P(n)) = f(P(n)) + f_el(n);
        end
    end
end


function FEMPlot(xVec,yVec,uVec,numRows)
    cla
    vx = uVec(1:2:end);
    vy = uVec(2:2:end);
    quiver(xVec,yVec,vx,vy,'LineWidth',1,'Color','k')
    hold on
    if length(numRows) > 1 % then tri els
        T = numRows;
        triplot(T,xVec,yVec,'Color','k')
    else
        numCols = length(xVec)/numRows;
        xMesh = reshape(xVec,[numCols,numRows])';
        yMesh = reshape(yVec,[numCols,numRows])';
        vxMat = reshape(vx,[numCols,numRows])';
        vyMat = reshape(vy,[numCols,numRows])';
        vMat = sqrt(vxMat.^2 + vyMat.^2);
        surf(xMesh,yMesh,zeros(size(vMat)),'FaceColor','none','EdgeColor','k')
    end
    view([0 90])
%     zlim([0 max(strainMat(:))+1])
%     caxis([0 max(strainMat(:))+1])
%     colorbar;
    daspect([1 1 1])
    drawnow
end