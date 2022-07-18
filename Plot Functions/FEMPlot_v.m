function FEMPlot_v(xVec,yVec,vx,vy,IEN)
%     cla
    % IEN is connectivity matrix
    quiver(xVec,yVec,vx,vy,'LineWidth',1,'Color','k')
    hold on
    if size(IEN,1) == 3
%         T = reshape(IEN,[3,length(IEN)/3])';
        T = IEN;
        triplot(T,xVec,yVec,'Color','k','LineWidth',0.1)
    else
%         numEl = size(IEN,1);
        numCols = IEN(4,1)-1; % for first element in the top left corner
        numRows = length(xVec)/numCols;
        xMesh = reshape(xVec,[numCols,numRows])';
        yMesh = reshape(yVec,[numCols,numRows])';
        vxMat = reshape(vx,[numCols,numRows])';
        vyMat = reshape(vy,[numCols,numRows])';
        vMat = sqrt(vxMat.^2 + vyMat.^2);
        surf(xMesh,yMesh,zeros(size(vMat)),'FaceColor','none','EdgeColor','k','LineWidth',0.1)
    end
    view([0 90])
%     zlim([0 max(strainMat(:))+1])
%     caxis([0 max(strainMat(:))+1])
%     colorbar;
    daspect([1 1 1])
    drawnow
end