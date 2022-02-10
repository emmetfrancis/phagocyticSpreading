function dtVal = dtCalc(rMesh,zMesh,v_rMesh,v_zMesh)
    numCols = sum(zMesh <= eps);
    numRows = length(rMesh)/numCols;
    numEl = (numCols-1) * (numRows-1);
    v_rInterp = scatteredInterpolant(rMesh,zMesh,v_rMesh);
    v_zInterp = scatteredInterpolant(rMesh,zMesh,v_zMesh);
    dtRefVals = zeros(1,numEl);
    for i = 1:numEl
        curRow = floor((i-1) / (numCols - 1)) + 1;
        curCol = mod(i-1, (numCols - 1)) + 1;
        startingNode = numCols * (curRow - 1) + curCol;
        endingNode = numCols * curRow + curCol;
        globalIdx = [startingNode,startingNode+1,endingNode+1,endingNode];
        rCur = rMesh(globalIdx);
        zCur = zMesh(globalIdx);
        curPolygon = polyshape(rMesh(globalIdx),zMesh(globalIdx));
        [rC,zC] = centroid(curPolygon);
        v_rCur = v_rInterp(rC,zC);
        v_zCur = v_zInterp(rC,zC);
        rDim = max(rCur) - min(rCur);
        zDim = max(zCur) - min(zCur);
        dtRefVals(i) = min([abs(rDim/v_rCur), abs(zDim/v_zCur)]);
    end
    dtVal = min(dtRefVals);
end