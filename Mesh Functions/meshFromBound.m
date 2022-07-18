function [vMeshStruct,pMeshStruct] = meshFromBound(r,z,nodesPerEl)%,topPlateLogic)
    % here, mesh is constructed via TFI with 
    % top: spaced out main cell body
    % right: concentrated points for adh calc
    % bottom: adherent cell
    % left: axis of symm
    % first, separate shape into top, right, bottom, left
    if size(r,1) > 1
        r = r';
        z = z';
    end
    sParam = zeros(1,length(r));
    for i = 2:length(sParam)
        sParam(i) = sParam(i-1) + sqrt((r(i)-r(i-1))^2 + (z(i)-z(i-1))^2);
    end
    surfIdx = find(z <= eps, 1, 'first');
    cornerIdx = find(sParam < sParam(surfIdx)/2,1,'last');
    numRows = surfIdx-cornerIdx + 1;
    numCols = cornerIdx;
    xiMesh = (0:1/(numCols-1):1)';
    etaMesh = (0:1/(numRows-1):1)';
    shapeCell{1} = [r(1:cornerIdx)',z(1:cornerIdx)'];
    shapeCell{2} = [r(surfIdx:-1:cornerIdx)',z(surfIdx:-1:cornerIdx)'];
    shapeCell{3} = [r(surfIdx)*xiMesh,zeros(size(xiMesh))];
    shapeCell{4} = [zeros(size(etaMesh)),z(1)*etaMesh];
    rNew = [r(1:surfIdx)'; shapeCell{3}(end-1:-1:1,1); shapeCell{4}(2:end-1,1)];
    zNew = [z(1:surfIdx)'; shapeCell{3}(end-1:-1:1,2); shapeCell{4}(2:end-1,2)];
    % make mesh for velocities
%     gSides_v = [ 0,0; 0,0; 0,0; 0,0 ]; %top,right,bottom,left
%     gSidesLogic_v = [0,0; 0,0; 1,1; 1,0]; % left is axis of symm here
%     hSides_v = [ 1,1; 1,1; 0,0; 0,0 ];
%     hSidesLogic_v = [ 1,1; 1,1; 0,0; 0,1 ];
%     gCell_v = {gSides_v,gSidesLogic_v};
%     hCell_v = {hSides_v,hSidesLogic_v};
    
    % make mesh for pressure
%     gSides_p = [ 0; 0; 0; 0 ]; %top,right,bottom,left
%     gSidesLogic_p = [0; 0; 0; 0]; % left is axis of symm here
%     hSides_p = [ 1; 1; 1; 1 ];
%     hSidesLogic_p = [ 1; 1; 1; 1 ];
%     gCell_p = {gSides_p,gSidesLogic_p};
%     hCell_p = {hSides_p,hSidesLogic_p};

    if nodesPerEl == 4
        [gCell_v,hCell_v] = setCellBC(rNew,zNew,'spread stick');
        [gCell_p,hCell_p] = setCellBC(rNew,zNew,'pressure');
        vMeshStruct = makeMeshNew(shapeCell,gCell_v,hCell_v);
        pMeshStruct = makeMeshNew(shapeCell,gCell_p,hCell_p);
    elseif nodesPerEl == 3
        [gCell_v,hCell_v] = setCellBC(r,z,'spread stick');
        [gCell_p,hCell_p] = setCellBC(r,z,'pressure');
        vMeshStruct = makeMeshDT(r',z',gCell_v,hCell_v);
        pMeshStruct = makeMeshDT(r',z',gCell_p,hCell_p);
    else
        error('invalid number of nodes per el')
    end
%     if topPlateLogic
%         offset = 5;
%         topIdx = find(z >= (max(z)-eps),1,'last') + offset;
%         cornerIdx = find(sParam < sParam(end)/4,1,'last');
% %         surfContactIdx = find(z <= eps, 1, 'first') - offset; 
%         if topIdx > cornerIdx %|| surfContactIdx < surfIdx
%             cornerIdx = topIdx;
% %             surfIdx = surfContactIdx;
%         end
%         surfIdx = length(sParam) - cornerIdx + 1;
%     end
%     if topPlateLogic
% %         gSides_v = [ 0,0; 0,0; 0,0; 0,0 ]; %top,right,bottom,left
% %         gSidesLogic_v = [0,0; 0,0; 0,1; 1,0]; % left is axis of symm here
% %         hSides_v = [ 0,1; 1,1; 0,0; 0,0 ];
% %         hSidesLogic_v = [ 1,1; 1,1; 1,0; 0,1 ];
%         numBound = 2*numCols + 2*numRows - 4;
%         gSides_v = zeros(2*numBound,1);
%         gSidesLogic_v = zeros(2*numBound,1);
%         hSides_v = zeros(2*numBound,1);
%         hSidesLogic_v = zeros(2*numBound,1);
%         topIdx =  find(z>=(max(z)-eps),1,'last');
%         surfIdx = find(z<=eps,1,'first');
%         leftIdx = numBound - numRows + 2;
%         gSidesLogic_v(2*surfIdx:2:2*leftIdx) = 1;
%         gSidesLogic_v(2*leftIdx-1:2:end-1) = 1;
% %         hSides_v(2:2:2*topIdx) = 1;
%         hSides_v(1:2*topIdx) = 1;
%         hSides_v(2*topIdx+1:2*(surfIdx-1)) = 1;
%         hSidesLogic_v(1:2*(surfIdx-1)) = 1;
%         hSidesLogic_v(2*surfIdx-1:2:2*leftIdx-1) = 1;
%         hSidesLogic_v(2*(leftIdx+1):2:end) = 1;
%     end
end