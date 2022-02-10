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
numRows = meshStruct.numRows;
numCols = meshStruct.numCols;
freeNodes = numRows * numCols - gNum;
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
numEl = (numRows-1)*(numCols-1);
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
volume = volCalc(x,y,axisymm);
% figure
if displayLogic
    FEMPlot(xVec,yVec,p0Vec,numRows)
end
% set(gcf, 'Position', [0 0 1000 1500])
% now calculate stiffness and mass matrices if not already computed
gVecFull = zeros(length(f),1);
for i = 1:numEl
    curIdx = (i-1) * 4 + 1;
    xCur = x(curIdx:curIdx+3);
    yCur = y(curIdx:curIdx+3);
    if ~QCalc
        K_el = elStiffness(xCur, yCur, axisymm);
        M_el = elMass(xCur, yCur, axisymm);
        K = KAssemble(P(curIdx:curIdx+3), K_el, K);
        M = MAssemble(P(curIdx:curIdx+3), M_el, M);
    end
%     g_fVec = K_el * gVec(curIdx:curIdx+3);
    gVecFull = fAssemble(P(curIdx:curIdx+3),gVec(curIdx:curIdx+3),gVecFull);
    fCell = {vx(curIdx:curIdx+3),vy(curIdx:curIdx+3),p0(curIdx:curIdx+3)/mu};
    f_elAlt = elForce(fCell, xCur, yCur, axisymm, true);
    f_el = elForce(fCell, xCur, yCur, axisymm, false);
%     fBC = elForceBC(hVec(curIdx:curIdx+3), hVal, xCur, yCur, axisymm);
%     %hVal=0 right now
    fBCDivergence = elForceBC(hVec(curIdx:curIdx+3),{vx(curIdx:curIdx+3),vy(curIdx:curIdx+3)}, xCur, yCur, axisymm);
    f_elAlt = f_elAlt + fBCDivergence;
    fAlt = fAssemble(P(curIdx:curIdx+3), f_elAlt, fAlt);
    f = fAssemble(P(curIdx:curIdx+3), f_el, f);
end
if ~QCalc
    Q = (1/mu)*M + (epsFactor/mu) * K;
    Q = sparse(Q);
    Qinv = inv(Q);
% else
%     Q = (1/mu)*M + epsilon * K;
%     Qinv = inv(Q);
end
p = Qinv * (f - (epsFactor/mu)*K*gVecFull);
% p = Q \ f;
% dimensionless
% pVec = p/(.025/8.5);
% xVec = xVec/8.5;
% yVec = yVec/8.5;
pVec = p;
if displayLogic
    FEMPlot(xVec,yVec,pVec,numRows)
end
QStruct{1} = K;
QStruct{2} = M;
QStruct{3} = Qinv;
end


function K_el = elStiffness(x, y, axisymm) 
    numNodes = 4;
    K_el = zeros(numNodes, numNodes);
    gaussPts = [-1/sqrt(3), 1/sqrt(3)];
    area = 0;
    for n = 1:2
        for m = 1:2
            curEps = gaussPts(n);
            curEta = gaussPts(m);
            N_eps = [ -(1 + curEta) / 4, (1 + curEta) / 4, (1 - curEta) / 4, -(1 - curEta) / 4 ]';
            N_eta = [ (1 - curEps) / 4, (1 + curEps) / 4, -(1 + curEps) / 4, -(1 - curEps) / 4 ]';
            x_eps = sum(N_eps .* x);
            x_eta = sum(N_eta .* x);
            y_eps = sum(N_eps .* y);
            y_eta = sum(N_eta .* y);
            j = x_eps * y_eta - x_eta * y_eps;
            jMat = (1/j) * [y_eta, -x_eta; -y_eps, x_eps];
            NMat = horzcat(N_eps, N_eta);
            NMat = NMat * jMat;
            if axisymm
                N_cur = [ 0.25 * (1 - curEps) * (1 + curEta), 0.25 * (1 + curEps) * (1 + curEta),...
                    0.25 * (1 + curEps) * (1 - curEta), 0.25 * (1 - curEps) * (1 - curEta) ]';
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
    numNodes = 4;
    fVec = zeros(numNodes,1);
    fVecAlt = zeros(numNodes,1);
    gaussPts = [-1/sqrt(3), 1/sqrt(3)];
    vx = fCell{1};
    vy = fCell{2};
    p0 = fCell{3};
    for n = 1:2
        for m = 1:2
            curEps = gaussPts(n);
            curEta = gaussPts(m);
            N_eps = [ -(1 + curEta) / 4, (1 + curEta) / 4, (1 - curEta) / 4, -(1 - curEta) / 4 ]';
            N_eta = [ (1 - curEps) / 4, (1 + curEps) / 4, -(1 + curEps) / 4, -(1 - curEps) / 4 ]';
            x_eps = sum(N_eps .* x);
            x_eta = sum(N_eta .* x);
            y_eps = sum(N_eps .* y);
            y_eta = sum(N_eta .* y);
            j = x_eps * y_eta - x_eta * y_eps;
            jMat = (1/j) * [y_eta, -x_eta; -y_eps, x_eps];
            NMat = horzcat(N_eps, N_eta);
            NMat = NMat * jMat;
            N_cur = [ 0.25 * (1 - curEps) * (1 + curEta), 0.25 * (1 + curEps) * (1 + curEta),...
                0.25 * (1 + curEps) * (1 - curEta), 0.25 * (1 - curEps) * (1 - curEta) ]';
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
    numNodes = 4;
    fVec = zeros(numNodes,1);
    epsVec = [-1, 1, 1, -1];
    etaVec = [1, 1, -1, -1];
    gaussPts = [-1/sqrt(3), 1/sqrt(3)];
    for n = 1:4
        if h(n) ~= 0
            if n == 4
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
    numNodes = 4;
    M_el = zeros(numNodes, numNodes);
    gaussPts = [-1/sqrt(3), 1/sqrt(3)];
    for n = 1:2
        for m = 1:2
            curEps = gaussPts(n);
            curEta = gaussPts(m);
            N_eps = [ -(1 + curEta) / 4, (1 + curEta) / 4, (1 - curEta) / 4, -(1 - curEta) / 4 ]';
            N_eta = [ (1 - curEps) / 4, (1 + curEps) / 4, -(1 + curEps) / 4, -(1 - curEps) / 4 ]';
            x_eps = sum(N_eps .* x);
            x_eta = sum(N_eta .* x);
            y_eps = sum(N_eps .* y);
            y_eta = sum(N_eta .* y);
            j = x_eps * y_eta - x_eta * y_eps;
            N_cur = [ 0.25 * (1 - curEps) * (1 + curEta), 0.25 * (1 + curEps) * (1 + curEta),...
                0.25 * (1 + curEps) * (1 - curEta), 0.25 * (1 - curEps) * (1 - curEta) ]';
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
    for n = 1:4
        if P(n) > 0
            for m = 1:4
                if P(m) > 0
                    K(P(n), P(m)) = K(P(n), P(m)) + K_el(n, m);
                end
            end
        end
    end
end


function M = MAssemble(P, M_el, M)
    for n = 1:4
        if P(n) > 0
            for m = 1:4
                if P(m) > 0
                    M(P(n), P(m)) = M(P(n), P(m)) + M_el(n, m);
                end
            end
        end
    end
end


function f = fAssemble(P, f_el, f)
    for n = 1:4
        if P(n) > 0
            f(P(n)) = f(P(n)) + f_el(n);
        end
    end
end

% function [P,IEN,x,y,gVec_global,hVec_global] = makeMesh(shapeCell, gCell, hCell)
% % rectangular domain split evenly o rectangular sections, 4 nodes each element
% % output is a (numRows x numCols + 2 x (numRows-2) + 2 x (numCols-2) + 3 x (numRows-2) x (numCols-2)) by 4 matrix, first column is P, second is x, third is y, fourth is gVec_global
% % EBC vector goes clockwise around the boundary   
%     %totNodes = numRows * numCols + 2 * (numRows - 2) + 2 * (numCols - 2) + 3 * (numRows - 2) * (numCols - 2);
%     xT = shapeCell{1};
%     xL = shapeCell{2};
%     xB = shapeCell{3};
%     xR = shapeCell{4};
%     numCols = length(xT);
%     numRows = length(xR);
%     numNodes = numRows * numCols;
%     [x,y,IEN] = TFI2DMesh(xT,xR,xB,xL);
%     x = x'; x = x(:);
%     y = y'; y = y(:);
%     IEN = IEN'; IEN = IEN(:);
%     numEl = (numCols - 1) * (numRows - 1);
%     gVec_lin = gCell{1};
%     gLogic_lin = gCell{2};
%     hVec_lin = hCell{1};
%     hLogic_lin = hCell{2};
%     gVec_global = zeros(4*numEl,1);
%     gLogic_global = zeros(4*numEl,1);
%     hVec_global = zeros(4*numEl,1);
%     for n = 1:numEl
%         startIdx = 4 * (n-1) + 1;
%         curRow = floor((n-1) / (numCols - 1)) + 1;
%         curCol = mod(n-1, (numCols - 1)) + 1;
%         startingNode = numCols * (curRow - 1) + curCol;
%         endingNode = numCols * curRow + curCol;
%         PIdx = startIdx:startIdx+3;
%         globalIdx = [startingNode,startingNode+1,endingNode+1,endingNode];
%         gIdx = zeros(1,4);
%         hIdx = zeros(1,4);
%         if curRow == 1
%             gIdx(1:2) = [curCol,curCol+1];
%             hIdx(1) = curCol;
%         elseif curRow == numRows-1
%             gIdx(3:4) = [(2*numCols + numRows - 1)-curCol-1,(2*numCols + numRows - 1)-curCol];
%             hIdx(3) = (2*numCols + numRows - 1)-curCol-1;
%         end
%         if curCol == 1
%             if curRow == 1 % to handle top left corner
%                 gIdx(4) = 2*numCols + 2*numRows - 2 - curRow-1;
%             else
%                 gIdx([4,1]) = [(2*numCols + 2*numRows - 2) - curRow-1,(2*numCols + 2*numRows - 2) - curRow];
%             end
%             hIdx(4) = (2*numCols + 2*numRows - 2) - curRow-1;
%         elseif curCol == numCols-1
%             gIdx(2:3) = [numCols + curRow - 1, numCols+curRow];
%             hIdx(2) = numCols + curRow - 1;
%         end
%         gVec_global(PIdx(gIdx ~= 0)) = gVec_lin(gIdx(gIdx ~= 0));
%         gLogic_global(PIdx(gIdx ~= 0)) = gLogic_lin(gIdx(gIdx~=0));
%         hVec_global(PIdx(hIdx ~= 0)) = hVec_lin(hIdx(hIdx ~= 0));
%     end
%     % assign equation numbers
%     curEqNum = 0;
%     P = zeros(4*numEl,1);
%     for n = 1:numNodes
%         elIdx = find(IEN == n);
%         if ~gLogic_global(elIdx(1))
%             curEqNum = curEqNum + 1;
%             P(elIdx) = curEqNum;
%         end
%     end
% end

function cSum = pTotCalc(x, y, P, c0, axisymm)
    numNodes = 4;
    cSum = 0;
    numEl = length(x) / numNodes;
    gaussPts = [ -1/sqrt(3), 1/sqrt(3) ];
    for i = 1:numEl
        curIdx = (i-1) * 4 + 1;
        xCur = x(curIdx:curIdx+3);
        yCur = y(curIdx:curIdx+3);
        cCur = c0(P(curIdx:curIdx+3));
        for n = 1:2
            for m = 1:2
                curEps = gaussPts(n);
                curEta = gaussPts(m);
                N_eps = [ -(1 + curEta) / 4, (1 + curEta) / 4, (1 - curEta) / 4, -(1 - curEta) / 4 ]';
                N_eta = [ (1 - curEps) / 4, (1 + curEps) / 4, -(1 + curEps) / 4, -(1 - curEps) / 4 ]';
                x_eps = sum(N_eps .* xCur);
                x_eta = sum(N_eta .* xCur);
                y_eps = sum(N_eps .* yCur);
                y_eta = sum(N_eta .* yCur);
                j = x_eps * y_eta - x_eta * y_eps;
                N_cur = [ 0.25 * (1 - curEps) * (1 + curEta), 0.25 * (1 + curEps) * (1 + curEta),...
                    0.25 * (1 + curEps) * (1 - curEta), 0.25 * (1 - curEps) * (1 - curEta) ]';
                if axisymm
                    xVal = sum(N_cur .* xCur);
                    curVal = j * 2*pi*xVal * sum(N_cur .* cCur);
                    cSum = cSum + curVal;
                else
                    curVal = j * sum(N_cur .* cCur);
                    cSum = cSum + curVal;
                end
            end
        end
    end
end

function vol = volCalc(x, y, axisymm)
    numNodes = 4;
    vol = 0;
    numEl = length(x) / numNodes;
    gaussPts = [ -1/sqrt(3), 1/sqrt(3) ];
    for i = 1:numEl
        curIdx = (i-1) * 4 + 1;
        xCur = x(curIdx:curIdx+3);
        yCur = y(curIdx:curIdx+3);
        for n = 1:2
            for m = 1:2
                curEps = gaussPts(n);
                curEta = gaussPts(m);
                N_eps = [ -(1 + curEta) / 4, (1 + curEta) / 4, (1 - curEta) / 4, -(1 - curEta) / 4 ]';
                N_eta = [ (1 - curEps) / 4, (1 + curEps) / 4, -(1 + curEps) / 4, -(1 - curEps) / 4 ]';
                x_eps = sum(N_eps .* xCur);
                x_eta = sum(N_eta .* xCur);
                y_eps = sum(N_eps .* yCur);
                y_eta = sum(N_eta .* yCur);
                j = x_eps * y_eta - x_eta * y_eps;
                N_cur = [ 0.25 * (1 - curEps) * (1 + curEta), 0.25 * (1 + curEps) * (1 + curEta),...
                    0.25 * (1 + curEps) * (1 - curEta), 0.25 * (1 - curEps) * (1 - curEta) ]';
                if axisymm
                    xVal = sum(N_cur .* xCur);
                    curVal = j * 2*pi*xVal;
                    vol = vol + curVal;
                else
                    curVal = j;
                    vol = vol + curVal;
                end
            end
        end
    end
end

function FEMPlot(xVec,yVec,cVec,numRows)
    numCols = length(cVec)/numRows;
    xMesh = reshape(xVec,[numCols,numRows]);
    yMesh = reshape(yVec,[numCols,numRows]);
    xMesh = xMesh';
    yMesh = yMesh';
    cMat = reshape(cVec,[numCols,numRows]);
    cMat = cMat';
    surf(xMesh,yMesh,cMat,cMat,'FaceColor','interp')
    view([0 90])
%     zlim([0 12])
%     caxis([0 12])
    colorbar;
    daspect([1 1 1])
    drawnow
end