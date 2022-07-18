function [pVec,QStruct,f] = FEMStokesPressureCalc(meshStruct,vxVec,vyVec,p0Vec,displayLogic,QStruct,epsFactor)
% initialize mesh
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
    allNodes = numRows*numCols;
    freeNodes = allNodes - gNum;
    nodesPerEl = 4;
    numEl = (numRows-1)*(numCols-1);
else % then triangles
    allNodes = max(IEN);
    freeNodes = allNodes - gNum;
    nodesPerEl = 3;
    numEl = meshStruct.numEl;
end

if isempty(QStruct)
    K = zeros(freeNodes, freeNodes);
    M = zeros(freeNodes, freeNodes);
    QCalc = false;
else
    K = QStruct{1};
    M = QStruct{2};
    Qinv = QStruct{3};
    QCalc = true;
end

% x = x*8.5;
% y = y*8.5;
xVec = zeros(freeNodes,1);
yVec = zeros(freeNodes,1);
for i = 1:length(IEN)
    xVec(IEN(i)) = x(i);
    yVec(IEN(i)) = y(i);
end
f = zeros(freeNodes,1);
fAlt = zeros(freeNodes,1);
% if isempty(vxVec) || isempty(vyVec) || isempty(p0Vec)
%     [coeffFile,path]= uigetfile('*.mat','Choose file for coefficient values');
%     coeffLoad = load(fullfile(path,coeffFile));
%     coeffVals = coeffLoad.coeffStored{2};
%     [vxVec,vyVec,p0Vec] = meshFlowCalc(xVec/8.5,yVec/8.5,coeffVals,'odd');
% %     vx = vx'; vy = vy'; p0 = p0';
% end
% vxVec = vxVec * .025/.2;
% vyVec = vyVec * .025/.2;
% p0Vec = p0Vec * .025/8.5;
vx = zeros(length(IEN),1);
vy = zeros(length(IEN),1);
p0 = zeros(length(IEN),1);
for i = 1:length(IEN)
%     p0Vec(IEN(i)) = p0(i);
    p0(i) = p0Vec(IEN(i));
    vx(i) = vxVec(IEN(i));
    vy(i) = vyVec(IEN(i));
end

mu = 1;
% epsilon = epsFactor * mean(abs(diff(x)))*mean(abs(diff(y))) / mu;
% epsilon = epsilon*.01;
hVal = 0;
totNodes = length(x);
% pTotInit = pTotCalc(x, y, P, p0,axisymm);
% figure
if displayLogic
    FEMPlot(xVec,yVec,p0Vec,numRows)
end
% set(gcf, 'Position', [0 0 1000 1500])
% now calculate stiffness and mass matrices if not already computed
gVecFull = zeros(length(f),1);
for i = 1:numEl
    curIdx = (i-1) * nodesPerEl + 1;
    xCur = x(curIdx:curIdx+(nodesPerEl-1));
    yCur = y(curIdx:curIdx+(nodesPerEl-1));
    if ~QCalc
        K_el = elStiffness(xCur, yCur, axisymm);
        M_el = elMass(xCur, yCur, axisymm);
        K = KAssemble(P(curIdx:curIdx+(nodesPerEl-1)), K_el, K);
        M = MAssemble(P(curIdx:curIdx+(nodesPerEl-1)), M_el, M);
    end
%     g_fVec = K_el * gVec(curIdx:curIdx+3);
    gVecFull = fAssemble(P(curIdx:curIdx+(nodesPerEl-1)),gVec(curIdx:curIdx+(nodesPerEl-1)),gVecFull);
    fCell = {vx(curIdx:curIdx+(nodesPerEl-1)),vy(curIdx:curIdx+(nodesPerEl-1)),p0(curIdx:curIdx+(nodesPerEl-1))/mu};
    f_elAlt = elForce(fCell, xCur, yCur, axisymm, true);
    f_el = elForce(fCell, xCur, yCur, axisymm, false);
%     fBC = elForceBC(hVec(curIdx:curIdx+3), hVal, xCur, yCur, axisymm);
%     %hVal=0 right now
    fBCDivergence = elForceBC(hVec(curIdx:curIdx+(nodesPerEl-1)),{vx(curIdx:curIdx+(nodesPerEl-1)),...
        vy(curIdx:curIdx+(nodesPerEl-1))}, xCur, yCur, axisymm);
    f_elAlt = f_elAlt + fBCDivergence;
    fAlt = fAssemble(P(curIdx:curIdx+(nodesPerEl-1)), f_elAlt, fAlt);
    f = fAssemble(P(curIdx:curIdx+(nodesPerEl-1)), f_el, f);
end
if ~QCalc
    Q = (1/mu)*M + (epsFactor/mu) * K;
%     Q = sparse(Q);
    Qinv = inv(Q);
    if any(isnan(Qinv(:))) || any(isinf(Qinv(:)))
        Qinv = pinv(Q);
    end
% else
%     Q = (1/mu)*M + (epsFactor/mu) * K;
%     Qinv = inv(Q);
end
p = Qinv * (f - (epsFactor/mu)*K*gVecFull);
% p = Q \ (f - (epsFactor/mu)*K*gVecFull);
% dimensionless
% pVec = p/(.025/8.5);
% xVec = xVec/8.5;
% yVec = yVec/8.5;
pVec = p;
if displayLogic
    switch nodesPerEl
        case 3
            FEMPlot(xVec,yVec,pVec,IEN)
        case 4
            FEMPlot(xVec,yVec,pVec,numRows)
    end
end
QStruct{1} = K;
QStruct{2} = M;
QStruct{3} = Qinv;
end


function K_el = elStiffness(x, y, axisymm) 
    numNodes = length(x);
    K_el = zeros(numNodes, numNodes);
    gaussPts = [-1/sqrt(3), 1/sqrt(3)];
    area = 0;
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
            j = x_eps * y_eta - x_eta * y_eps;
            jMat = (1/j) * [y_eta, -x_eta; -y_eps, x_eps];
            NMat = horzcat(N_eps, N_eta);
            NMat = NMat * jMat;
            if axisymm
                xCur = sum(N_cur .* x);
                K_el = K_el + j * xCur * (NMat * NMat');
            else
                K_el = K_el + j * (NMat * NMat');
            end
            area = area + j;
        end
    end
%     epsCur = (max(x)-min(x))*(max(y)-min(y))/4;
    epsCur = area/pi;
    K_el = epsCur * K_el;
end

function fVec = elForce(fCell, x, y, axisymm, altCalc)
    numNodes = length(x);
    fVec = zeros(numNodes,1);
    fVecAlt = zeros(numNodes,1);
    gaussPts = [-1/sqrt(3), 1/sqrt(3)];
    vx = fCell{1};
    vy = fCell{2};
    p0 = fCell{3};
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
            j = x_eps * y_eta - x_eta * y_eps;
            jMat = (1/j) * [y_eta, -x_eta; -y_eps, x_eps];
            NMat = horzcat(N_eps, N_eta);
            NMat = NMat * jMat;
            pCur = sum(N_cur .* p0);
            vxCur = sum(N_cur .* vx);
            vyCur = sum(N_cur .* vy);
            if axisymm
                xCur = sum(N_cur .* x);
                fVecAlt = fVecAlt + j*xCur*(pCur*N_cur + vxCur*NMat(:,1) + vyCur*NMat(:,2));
                NMat(:,1) = NMat(:,1) + N_cur/xCur;
                vxDer = sum(NMat(:,1) .* vx);
                vyDer = sum(NMat(:,2) .* vy);
                f = vxDer + vyDer - pCur;
                fVec = fVec - f*j*xCur*N_cur;
            else
%                 vxDer = sum(NMat(:,1) .* vx);
%                 vyDer = sum(NMat(:,2) .* vy);
%                 f = vxDer + vyDer - pCur;
%                 fVec = fVec - f*j*N_cur;
                fVec = fVec + j*(pCur*N_cur + vxCur*NMat(:,1) + vyCur*NMat(:,2));
            end
        end
    end
%     fVecAlt(1:2) = -fVecAlt(1:2);
    if altCalc
        fVec = fVecAlt;
    end
    % finally subtract contribution from EBC
%     fVec = fVec - g_fVec;
end


function fVec = elForceBC(h,hVal,x,y,axisymm)
% for diffusion problem, compute boundary flux based on previous c
% need to do a 1d integral along each side with a boundary condition
    numNodes = length(x);
    fVec = zeros(numNodes,1);
    epsVec = [-1, 1, 1, -1];
    etaVec = [1, 1, -1, -1];
    gaussPts = [-1/sqrt(3), 1/sqrt(3)];
    for n = 1:numNodes
        if h(n) ~= 0
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
                N_cur = [ 0.25 * (1 - curEps) * (1 + curEta), 0.25 * (1 + curEps) * (1 + curEta),...
                        0.25 * (1 + curEps) * (1 - curEta), 0.25 * (1 - curEps) * (1 - curEta) ]';
                    if numNodes == 3
                        N_cur = [N_cur(1:2); N_cur(3)+N_cur(4)];
                    end
%                 N_eps = [ -(1 + curEta) / 4, (1 + curEta) / 4, (1 - curEta) / 4, -(1 - curEta) / 4 ]';
%                 N_eta = [ (1 - curEps) / 4, (1 + curEps) / 4, -(1 + curEps) / 4, -(1 - curEps) / 4 ]';
%                 c_eps = sum(N_eps .* c0);
%                 c_eta = sum(N_eta .* c0);
                dx = x(nextIdx)-x(n);
                dy = y(nextIdx)-y(n);
                dL = sqrt(dx^2 + dy^2);
                j1D = dL / 2;
                eps1D = gaussPts(m);
                N1DVal = 0.5*[1-eps1D;1+eps1D];
                if iscell(hVal)
                    n_x = -dy / dL;
                    n_y = dx / dL;
                    vxVal = sum(N_cur .* hVal{1});
                    vyVal = sum(N_cur .* hVal{2});
                    hCur = n_x * vxVal + n_y * vyVal;
                else
                    hCur = hVal;
                end
%                 x_eps = sum(N_eps .* x);
%                 x_eta = sum(N_eta .* x);
%                 y_eps = sum(N_eps .* y);
%                 y_eta = sum(N_eta .* y);
%                 j2D = x_eps * y_eta - x_eta * y_eps;
%                 jMat = (1/j2D) * [y_eta, -x_eta; -y_eps, x_eps];
%                 cDer = [c_eps, c_eta] * jMat;
%                 c_x = cDer(1);
%                 c_y = cDer(2);
%                 hCur = c_x * n_x + c_y * n_y;
%                 if hCur > 1
%                     fprintf('too large')
%                 end
                if axisymm
                    xCur = sum(N_cur .* x);
                    fVec(n) = fVec(n) - j1D * xCur * hCur * N1DVal(1);
                    fVec(nextIdx) = fVec(nextIdx) - j1D * xCur * hCur * N1DVal(2);
                else
                    fVec(n) = fVec(n) - j1D * hCur * N1DVal(1);
                    fVec(nextIdx) = fVec(nextIdx) - j1D * hCur * N1DVal(2);
                end
            end
        end
    end
end


function M_el = elMass(x, y, axisymm)
    numNodes = length(x);
    M_el = zeros(numNodes, numNodes);
    gaussPts = [-1/sqrt(3), 1/sqrt(3)];
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
            j = x_eps * y_eta - x_eta * y_eps;
            if axisymm
                xCur = sum(N_cur .* x);
                M_el = M_el + j * xCur * (N_cur * N_cur');
            else
                M_el = M_el + j * (N_cur * N_cur');
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


function M = MAssemble(P, M_el, M)
    for n = 1:length(P)
        if P(n) > 0
            for m = 1:length(P)
                if P(m) > 0
                    M(P(n), P(m)) = M(P(n), P(m)) + M_el(n, m);
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

function FEMPlot(xVec,yVec,cVec,numRows)
    if length(numRows) > 1 % then numRows is actually triangle connectivity
        T = numRows;
        trisurf(T,xVec,yVec,cVec,cVec,'FaceColor','interp')
    else
        numCols = length(cVec)/numRows;
        xMesh = reshape(xVec,[numCols,numRows]);
        yMesh = reshape(yVec,[numCols,numRows]);
        cMat = reshape(cVec,[numCols,numRows]);
        xMesh = xMesh';
        yMesh = yMesh';
        cMat = cMat';
        surf(xMesh,yMesh,cMat,cMat,'FaceColor','interp')
    end
    view([0 90])
%     zlim([0 12])
%     caxis([0 12])
    colorbar;
    daspect([1 1 1])
    drawnow
end