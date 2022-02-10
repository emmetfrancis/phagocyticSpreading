function [xMesh,yMesh,IEN] = TFI2DMesh(xT,xR,xB,xL)
% Cuts shape into quadrilaterals using TFI (transfinite interpolation), 
% store x and y in matrices (each row is an element), store connectivity
% INPUTS:
% xT: Nx by 2 matrix of x and y values along the "top" (xi = 0->1, eta = 1)
% xR: Ny by 2 matrix of x and y values along the "right" (xi = 1, eta = 0->1)
% xB: Nx by 2 matrix of x and y values along the "bottom" (xi = 0->1, eta = 0)
% xL: Ny by 2 matrix of x and y values along the "left" (xi = 0, eta = 0->1)
% xi is length Nx and eta is length Ny, Nx and Ny set fineness of mesh (numEl = (Nx-1)*(Ny-1))
% OUTPUTS:
% xMesh: (Nx-1)*(Ny-1) by 4 matrix of x values (every four defines an element)
% yMesh: (Nx-1)*(Ny-1) by 4 matrix of y values (every four defines an element)
% IEN: local to global node numbers ((Nx-1)*(Ny-1) by 4) (connectivity)
[Nx,whichNx] = max([size(xT,1),size(xB,1)]);
xi = 0:1/(Nx-1):1;
[Ny,whichNy] = max([size(xR,1),size(xL,1)]);
eta = 0:1/(Ny-1):1;
if size(xT,1) ~= size(xB,1)
    switch whichNx
        case 1
            otherXi = 0:1/(size(xB,1)-1):1;
            xB = interp1q(otherXi,xB,xi);
        case 2
            otherXi = 0:1/(size(xT,1)-1):1;
            xT = interp1q(otherXi,xT,xi);
    end
end
if size(xR,1) ~= size(xL,1)
    switch whichNy
        case 1
            otherEta = 0:1/(size(xL,1)-1):1;
            xL = interp1q(otherEta,xL,eta);
        case 2
            otherEta = 0:1/(size(xR,1)-1):1;
            xR = interp1q(otherEta,xR,eta);
    end
end
% format mesh in logical coordinates
xiMesh = zeros((Nx-1)*(Ny-1),4);
etaMesh = zeros((Nx-1)*(Ny-1),4);
dXi = 1/(Nx-1);
dEta = 1/(Ny-1);
IEN = zeros((Nx-1)*(Ny-1),4);
for i = 1:(Nx-1)*(Ny-1)
    curRow = ceil(i/(Nx-1));
    curCol = mod((i-1),Nx-1) + 1;
    xi1 = (curCol-1)/(Nx-1);
    xiMesh(i,:) =  [xi1, xi1 + dXi, xi1 + dXi, xi1];
    eta1 = 1-(curRow-1)/(Ny-1);
    etaMesh(i,:) = [eta1, eta1, eta1-dEta, eta1-dEta];
    IEN(i,:) = [Nx*(curRow-1)+curCol, Nx*(curRow-1)+curCol+1,...
        Nx*curRow+curCol+1, Nx*curRow+curCol];
end
xIdxMesh = round(xiMesh*(Nx-1)+1);
yIdxMesh = round(etaMesh*(Ny-1)+1);
% TFI calc
xMesh = (1-etaMesh) .* xB(xIdxMesh) + etaMesh .* xT(xIdxMesh) +...
    (1-xiMesh) .* xL(yIdxMesh) + xiMesh .* xR(yIdxMesh) -...
    (xiMesh .* etaMesh .* xT(end,1) + xiMesh .* (1-etaMesh) .* xB(end,1) +...
    etaMesh .* (1-xiMesh) .* xT(1,1) + (1-xiMesh) .* (1-etaMesh) .* xB(1,1));
yMesh = (1-etaMesh) .* xB(Nx+xIdxMesh) + etaMesh .* xT(Nx+xIdxMesh) +...
    (1-xiMesh) .* xL(Ny+yIdxMesh) + xiMesh .* xR(Ny+yIdxMesh) -...
    (xiMesh .* etaMesh .* xT(end,2) + xiMesh .* (1-etaMesh) .* xB(end,2) +...
    etaMesh .* (1-xiMesh) .* xT(1,2) + (1-xiMesh) .* (1-etaMesh) .* xB(1,2));
% node numbering
% curNum = 0;
% IEN = zeros((Nx-1)*(Ny-1),4);
% xLin = ;
% for i = 1:size(IEN,1)
%     % local numbers of row are 1 to 4
%     curIdx = i*4-3:i*4;
%     pastIdx = 1:i*4-4;
%     for j = 1:4
%         curLogic = (xLin(pastIdx) == xLin(curIdx(j)) & yLin(pastIdx) == yLin(curIdx(j)));
%         if any(curLogic)
%             linIdx = find(curLogic,1,'first');
%             IENTrans = IEN';
%             IEN(i,j) = IENTrans(linIdx);
%         else
%             IEN(i,j) = curNum + 1;
%             curNum = curNum + 1; % number of global nodes so far
%         end
%     end
% end
end