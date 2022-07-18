function FEMPlot_p(xVec,yVec,pVec,IEN)
    % IEN is actually connectivity matrix
    if size(IEN,1) == 3
%         T = reshape(IEN,[3,length(IEN)/3])';
        T = IEN;
        trisurf(T,xVec,yVec,pVec,pVec,'FaceColor','interp')
    else
%         numEl = size(IEN,1);
        numCols = IEN(4,1)-1; % for first element in the top left corner
        numRows = length(xVec)/numCols;
        xMesh = reshape(xVec,[numCols,numRows]);
        yMesh = reshape(yVec,[numCols,numRows]);
        xMesh = xMesh';
        yMesh = yMesh';
        pMat = reshape(pVec,[numCols,numRows]);
        pMat = pMat';
        surf(xMesh,yMesh,pMat,pMat,'FaceColor','interp')
    end
    view([0 90])
%     zlim([-50 20])
%     caxis([-50 20])
    colorbar;
    daspect([1 1 1])
    drawnow
end