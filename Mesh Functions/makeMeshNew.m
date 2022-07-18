function meshStruct = makeMeshNew(shapeCell, gCell, hCell)
% rectangular domain split evenly o rectangular sections, 4 nodes each element
% output is a (numRows x numCols + 2 x (numRows-2) + 2 x (numCols-2) + 3 x (numRows-2) x (numCols-2)) by 4 matrix, first column is P, second is x, third is y, fourth is gVec_global
% EBC vector goes clockwise around the boundary   
    %totNodes = numRows * numCols + 2 * (numRows - 2) + 2 * (numCols - 2) + 3 * (numRows - 2) * (numCols - 2);
    axisymm = true;
    xT = shapeCell{1};
    xR = shapeCell{2};
    xB = shapeCell{3};
    xL = shapeCell{4};
    numCols = length(xT);
    numRows = length(xR);
    boundLength = 2*numCols + 2*(numRows-2);
    if length(gCell{1}) == 2*boundLength || length(gCell{1}) == boundLength
        gVec_lin = gCell{1};
        gLogic_lin = gCell{2};
        gDOF = length(gVec_lin) / boundLength;
        hVec_lin = hCell{1};
%         hLogic_lin = hCell{2};
        hDOF = length(hVec_lin) / boundLength;
    else
        gSides = gCell{1}; %top,right,bottom,left
        gSidesLogic = gCell{2};
        hSides = hCell{1};
        hSidesLogic = hCell{2};
        gDOF = size(gSides,2);
        hDOF = size(hSides,2);
        gVec_lin = zeros(gDOF*boundLength,1);
        gLogic_lin = zeros(gDOF*boundLength,1);
        hVec_lin = zeros(hDOF*boundLength,1);
        hLogic_lin = zeros(hDOF*boundLength,1);
        % construct the boundary condition vectors
        for n = 1:boundLength
            % on the corners, g will be averaged from the side values
            if n < numCols
                if n == 1
                    gVec_lin(gDOF*(n-1)+1:gDOF*n) = mean([gSides(1,:); gSides(4,:)]);
                    gLogic_lin(gDOF*(n-1)+1:gDOF*n) = gSidesLogic(1,:) | gSidesLogic(4,:);
                else
                    gVec_lin(gDOF*(n-1)+1:gDOF*n) = gSides(1,:);
                    gLogic_lin(gDOF*(n-1)+1:gDOF*n) = gSidesLogic(1,:);
                end
                hVec_lin(hDOF*(n-1)+1:hDOF*n) = hSides(1,:);
                hLogic_lin(hDOF*(n-1)+1:hDOF*n) = hSidesLogic(1,:);
            elseif n < numCols+numRows-1
                if n == numCols
                    gVec_lin(gDOF*(n-1)+1:gDOF*n) = mean([gSides(1,:); gSides(2,:)]);
                    gLogic_lin(gDOF*(n-1)+1:gDOF*n) = gSidesLogic(1,:) | gSidesLogic(2,:);
                else
                    gVec_lin(gDOF*(n-1)+1:gDOF*n) = gSides(2,:);
                    gLogic_lin(gDOF*(n-1)+1:gDOF*n) = gSidesLogic(2,:);
                end
                hVec_lin(hDOF*(n-1)+1:hDOF*n) = hSides(2,:);
                hLogic_lin(hDOF*(n-1)+1:hDOF*n) = hSidesLogic(2,:);
            elseif n<(2*numCols+numRows-2)
                if n == numCols+numRows-1
                    gVec_lin(gDOF*(n-1)+1:gDOF*n) = mean([gSides(2,:); gSides(3,:)]);
                    gLogic_lin(gDOF*(n-1)+1:gDOF*n) = gSidesLogic(2,:) | gSidesLogic(3,:);
                else
                    gVec_lin(gDOF*(n-1)+1:gDOF*n) = gSides(3,:);
                    gLogic_lin(gDOF*(n-1)+1:gDOF*n) = gSidesLogic(3,:);
                end
                hVec_lin(hDOF*(n-1)+1:hDOF*n) = hSides(3,:);
                hLogic_lin(hDOF*(n-1)+1:hDOF*n) = hSidesLogic(3,:);
            else
                if n == (2*numCols+numRows-2)
                    gVec_lin(gDOF*(n-1)+1:gDOF*n) = mean([gSides(3,:); gSides(4,:)]);
                    gLogic_lin(gDOF*(n-1)+1:gDOF*n) = gSidesLogic(3,:) | gSidesLogic(4,:);
                else
                    gVec_lin(gDOF*(n-1)+1:gDOF*n) = gSides(4,:);
                    gLogic_lin(gDOF*(n-1)+1:gDOF*n) = gSidesLogic(4,:);
                end
                hVec_lin(hDOF*(n-1)+1:hDOF*n) = hSides(4,:);
                hLogic_lin(hDOF*(n-1)+1:hDOF*n) = hSidesLogic(4,:);
            end
        end
        % gvec is zeros, hvec is 0 initially but needs to be computed after each time step
    end
    gNum = sum(gLogic_lin);
    numNodes = numRows * numCols;
    numEl = (numCols - 1) * (numRows - 1);
    [x,y,IEN] = TFI2DMesh(xT,xR,xB,xL);
    x = x'; x = x(:);
    y = y'; y = y(:);
    IEN = IEN'; IEN = IEN(:);
    nDOF = length(gVec_lin)/boundLength; % dof per node
    if nDOF ~= 1 && nDOF ~= 2
        error('mismatch in dimensions for mesh vs. gVec')
    end
    hDOF = length(hVec_lin)/boundLength;
    if hDOF ~= 1 && hDOF ~= 2
        error('mismatch in dimensions for mesh vs. hVec')
    end
    
    gVec_global = zeros(4*nDOF*numEl,1);
    gLogic_global = zeros(4*nDOF*numEl,1);
    hVec_global = zeros(4*hDOF*numEl,1);
    
    for n = 1:numEl
        PStartIdx = 4*nDOF * (n-1) + 1;
        hStartIdx = 4*hDOF * (n-1) + 1;
        xStartIdx = 4 * (n-1) + 1;
        curRow = floor((n-1) / (numCols - 1)) + 1;
        curCol = mod(n-1, (numCols - 1)) + 1;
        startingNode = numCols * (curRow - 1) + curCol;
        endingNode = numCols * curRow + curCol;
        PIdx = PStartIdx:PStartIdx+4*nDOF-1;
        hGlobalIdx = hStartIdx:hStartIdx+4*hDOF-1;
        xIdx = xStartIdx:xStartIdx+3;
        globalIdx = [startingNode,startingNode+1,endingNode+1,endingNode];
        gIdx = zeros(1,4*nDOF);
        hIdx = zeros(1,4*hDOF);
        if curRow == 1
            gIdx(1:2*nDOF) = nDOF*(curCol-1)+1:nDOF*(curCol+1);
            hIdx(1:hDOF) = hDOF*(curCol-1)+1:hDOF*(curCol);
            xEdge(curCol) = x(1);
        elseif curRow == numRows-1
            edgeIdx = (2*numCols + numRows - 1)-curCol-1;
            gIdx(2*nDOF+1:4*nDOF) = nDOF*(edgeIdx-1)+1:nDOF*(edgeIdx+1); 
            hIdx(2*hDOF+1:3*hDOF) = hDOF*(edgeIdx-1)+1:hDOF*edgeIdx;
        end
        if curCol == 1
            edgeIdx = 2*numCols + 2*numRows - 2 - curRow-1;
            hIdx(3*hDOF+1:4*hDOF) = hDOF*(edgeIdx-1)+1:hDOF*edgeIdx;
            if curRow == 1 % to handle top left corner
                gIdx(3*nDOF+1:4*nDOF) = nDOF*(edgeIdx-1)+1:nDOF*edgeIdx;
            else 
                gIdx([3*nDOF+1:4*nDOF,1:nDOF]) = nDOF*(edgeIdx-1)+1:nDOF*(edgeIdx+1);
            end
        elseif curCol == numCols-1
            edgeIdx = numCols + curRow - 1;
            gIdx(nDOF+1:3*nDOF) = nDOF*(edgeIdx-1)+1:nDOF*(edgeIdx+1);
            hIdx(hDOF+1:2*hDOF) = hDOF*(edgeIdx-1)+1:hDOF*edgeIdx;
        end
        gVec_global(PIdx(gIdx ~= 0)) = gVec_lin(gIdx(gIdx ~= 0));
        gLogic_global(PIdx(gIdx ~= 0)) = gLogic_lin(gIdx(gIdx~=0));
        hVec_global(hGlobalIdx(hIdx ~= 0)) = hVec_lin(hIdx(hIdx ~= 0));
    end
    
    % assign equation numbers
    curEqNum = 0;
    P = zeros(4*nDOF*numEl,1);
    for n = 1:numNodes
        elIdx = find(IEN == n);
        for i = 1:nDOF
            if ~gLogic_global(nDOF*(elIdx(1)-1)+i)
                curEqNum = curEqNum + 1;
                P(nDOF*(elIdx-1)+i) = curEqNum;
            end
        end
    end
    
    xVec = zeros(numRows*numCols,1);
    yVec = zeros(numRows*numCols,1);
    for i = 1:length(IEN)
        xVec(IEN(i)) = x(i);
        yVec(IEN(i)) = y(i);
    end
    
    % define boundary nodes
    topIndices = 1:numCols-1;
    rightIndices = numCols*(1:numRows-1);
    bottomIndices = numRows*numCols:-1:((numRows-1)*numCols+2);
    leftIndices = numCols*(numRows-1:-1:1)+1;
    meshStruct.boundIndices = [topIndices,rightIndices,bottomIndices,leftIndices(1)];
    
    %output into struct
    meshStruct.nodesPerEl = 4;
    meshStruct.P = P;
    meshStruct.IEN = IEN;
    meshStruct.x = x;
    meshStruct.y = y;
    meshStruct.gVec = gVec_global;
    meshStruct.hVec = hVec_global;
    meshStruct.gNum = gNum;
    meshStruct.axisymm = axisymm;
    meshStruct.numRows = numRows;
    meshStruct.numCols = numCols;
    meshStruct.xVec = xVec;
    meshStruct.yVec = yVec;
    
end