function [xBound,yBound] = meshToBound(xMesh,yMesh)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    numCols = sum(yMesh<eps);
    numRows = length(yMesh)/numCols;
    topIndices = 1:numCols-1;
    rightIndices = numCols*(1:numRows-1);
    bottomIndices = numRows*numCols:-1:((numRows-1)*numCols+2);
    leftIndices = numCols*(numRows-1:-1:1)+1;
    boundIndices = [topIndices,rightIndices,bottomIndices,leftIndices(1)];
    xBound = xMesh(boundIndices)';
    yBound = yMesh(boundIndices)';
end

