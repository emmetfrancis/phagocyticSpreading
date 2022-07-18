function meshStruct = makeMeshDT(r,z,gCell,hCell)
% [r,z]: boundary coordinates
% gLogic: length is 2*length(r), sets whether each node is fixed by EBC
gVec_lin = gCell{1};
gLogic_lin = gCell{2};
hVec_lin = hCell{1};
boundLength = length(r);
boundIndices = 1:length(r);
gDOF = length(gVec_lin) / boundLength;
% hLogic_lin = hCell{2}; % hLogic can be inferred from gLogic
hDOF = length(hVec_lin) / boundLength;
if gDOF ~= 1 && gDOF ~= 2
    error('mismatch in dimensions for mesh vs. gVec')
end
if hDOF ~= 1 && hDOF ~= 2
    error('mismatch in dimensions for mesh vs. hVec')
end

numLeftNodes = round((z(1)-z(end))/.02);
leftInt = (z(1)-z(end))/numLeftNodes;
r = [r; zeros(numLeftNodes-1,1)];
z = [z; (leftInt:leftInt:z(1)-leftInt)'];
boundVec = [r,z];
gVec_lin = [gVec_lin; zeros(gDOF*(numLeftNodes-1),1)];
gLogic_lin = [gLogic_lin; zeros(gDOF*(numLeftNodes-1),1)];
if gDOF == 2
    gLogic_lin(end-1:-2:end-2*(numLeftNodes-1)-1) = 1;
end
hVec_lin = [hVec_lin; ones(hDOF*(numLeftNodes-1),1)];
if hDOF == 2
    hVec_lin(end-1:-2:end-2*(numLeftNodes-1)-1) = 0;
end
gBoundIdx = 1:gDOF*length(r);
hBoundIdx = 1:hDOF*length(r);

sParam = zeros(size(boundVec,1),1);
for i = 2:size(boundVec,1)
    sParam(i) = sParam(i-1) + norm(boundVec(i,:)-boundVec(i-1,:));
end
% cx = (1/sParam(end)) * trapz(sParam,boundVec(:,1)); % "contour centroid"
% cy = (1/sParam(end)) * trapz(sParam,boundVec(:,2));
% % 2 boundary layers
startSurf = find(z<=eps,1,'first');
rBody = r(1:startSurf);
zBody = z(1:startSurf);
[~,~,phiValsBody] = polyCurv([],rBody,zBody,2,7);
bottomLeftIdx = find(r(2:end)<=eps,1,'first') + 1;
phiVals = [-pi/4; phiValsBody(2:end)'; pi*ones(bottomLeftIdx-startSurf-1,1);...
    5*pi/4; 3*(pi/2)*ones(length(r)-bottomLeftIdx,1)];
innerBoundVec1 = [r-.005*sin(phiVals), z-.005*cos(phiVals)];
spaceInt = sParam(end)/(round(sParam(end)/.02));
sParamSpaced = 0:spaceInt:sParam(end);
rSpaced = interp1(sParam,r,sParamSpaced)';
zSpaced = interp1(sParam,z,sParamSpaced)';
phiSpaced = interp1(sParam,phiVals,sParamSpaced)';
% innerBoundVec2 = [rSpaced-.02*sin(phiSpaced), zSpaced-.02*cos(phiSpaced)];
innerBoundVec2 = [rSpaced-.01*sin(phiSpaced), zSpaced-.01*cos(phiSpaced)];

% innerBoundVec1 = .92*[boundVec(:,1)-cx, boundVec(:,2)-cy] + [cx,cy];
% innerBoundVec2 = .85*[boundVec(:,1)-cx, boundVec(:,2)-cy] + [cx,cy];
% adaptively spaced, refined mesh towards surface
meshSpacing = .02; % for main body, where not refined
xMin = min(innerBoundVec2(:,1));
xMax = max(innerBoundVec2(:,1));
xNumInt = round((xMax-xMin)/meshSpacing);
xInt = (xMax-xMin)/xNumInt;
xMesh = xMin+xInt:xInt:xMax-xInt;
% yMin = min(boundVec(:,2));
% if max(boundVec(:,2))<2/8.5
%     mainMesh = [];
% else
%     yMin = min([2/8.5, max(boundVec(:,2))]);
    yMin = min(innerBoundVec2(:,2));
    yMax = max(innerBoundVec2(:,2));
    yNumInt = round((yMax-yMin)/meshSpacing);
    yInt = (yMax-yMin)/yNumInt;
    yMesh = yMin:yInt:yMax-yInt;
    [xMesh,yMesh] = meshgrid(xMesh,yMesh);
    for i = 2:2:size(xMesh,1)
        xMesh(i,:) = xMesh(i,:) + meshSpacing/2;
    end
    xMesh = xMesh(:);
    yMesh = yMesh(:);
    in = inpolygon(xMesh,yMesh,boundVec(:,1),boundVec(:,2));
    xMesh = xMesh(in);
    yMesh = yMesh(in);
    mainMesh = [xMesh,yMesh];
% end

% startSurf = find(z < eps, 1,'first');
% leftIdx = find(r(2:end) < eps,1,'first') + 1;
% xRefMin = boundVec(startSurf,1);
% xRefMax = max(boundVec(:,1));
% xNumRefInt = round((xRefMax-xRefMin)/(meshSpacing/2));
% xRefInt = (xRefMax-xRefMin)/xNumRefInt;
% xRefMesh = [boundVec(leftIdx-1:-1:startSurf,1)',xRefMin+xRefInt:xRefInt:xRefMax-xRefInt];
% yRefStart = find(z < 2/8.5,1,'first');
% yRefEnd = startSurf - 1;
% yRefMesh = z(yRefEnd:-1:yRefStart);
% [xRefMesh,yRefMesh] = meshgrid(xRefMesh,yRefMesh);
% for i = 2:2:size(xRefMesh,1)
%     xRefMesh(i,:) = xRefMesh(i,:) + meshSpacing/4;
% end
% xRefMesh = xRefMesh(:);
% yRefMesh = yRefMesh(:);
% in = inpolygon(xRefMesh,yRefMesh,boundVec(:,1),boundVec(:,2));
% xRefMesh = xRefMesh(in);
% yRefMesh = yRefMesh(in);
% refinedMesh = [xRefMesh,yRefMesh];

% allVertices = vertcat(boundVec, innerBoundVec1, innerBoundVec2, mainMesh);
allVertices = vertcat(boundVec, innerBoundVec2, mainMesh);
% allVertices = vertcat(boundVec, mainMesh, refinedMesh);
edgeID = zeros(length(boundVec),2);
edgeID(:,1) = 1:length(boundVec);
edgeID(1:end-1,2) = 2:length(boundVec);
edgeID(end,2) = 1;
DT = delaunayTriangulation(allVertices,edgeID);
IEN = DT.ConnectivityList;
% newIEN = [1,size(boundVec,1),2]; % initialize with top el
newIEN = [];
for i = 1:size(IEN,1)
   if ~all(IEN(i,:)<=size(boundVec,1))
        newIEN = [newIEN; IEN(i,:)]; % prohibit triangles outside cell body
   end
end
IEN = newIEN;
xVec = allVertices(:,1);
yVec = allVertices(:,2);
x = zeros(size(IEN));
y = zeros(size(IEN));
for i = 1:size(IEN,1) 
    IEN(i,:) = IEN(i,end:-1:1); % redefine counterclockwise
    x(i,:) = xVec(IEN(i,:));
    y(i,:) = yVec(IEN(i,:));
end

gNum = sum(gLogic_lin);
numNodes = length(xVec);
numEl = size(IEN,1);
x = x'; x = x(:);
y = y'; y = y(:);
IEN = IEN'; IEN = IEN(:);

gVecMesh = zeros(gDOF*size(xVec,1),1);
gVecMesh(gBoundIdx) = gVec_lin;
gLogicMesh = zeros(gDOF*length(xVec),1);
gLogicMesh(gBoundIdx) = gLogic_lin;
hVecMesh = zeros(hDOF*length(xVec),1);
hVecMesh(hBoundIdx) = hVec_lin;

gVec_global = zeros(3*gDOF*numEl,1);
gLogic_global = zeros(3*gDOF*numEl,1);
hVec_global = zeros(3*hDOF*numEl,1);

% assign equation numbers
curEqNum = 0;
P = zeros(3*gDOF*numEl,1);
for n = 1:numNodes
    elIdx = find(IEN == n);
    gVecCur = gVecMesh(gDOF*(n-1)+1:gDOF*n);
    gLogicCur = gLogicMesh(gDOF*(n-1)+1:gDOF*n);
    hVecCur = hVecMesh(hDOF*(n-1)+1:hDOF*n);
    for i = 1:gDOF
        gIdx = gDOF*(elIdx-1)+i;
        gVec_global(gIdx) = gVecCur(i);
        gLogic_global(gIdx) = gLogicCur(i);
        if ~gLogicCur(i)
            curEqNum = curEqNum + 1;
            P(gDOF*(elIdx-1)+i) = curEqNum;
        end
    end  
end
for n = 1:numNodes
    elIdx = find(IEN == n);
    hVecCur = hVecMesh(hDOF*(n-1)+1:hDOF*n);
    for j = 1:length(elIdx)
        localIdx = mod(elIdx(j),3);
        if localIdx == 0
            nextIdxGlobal = elIdx(j)-2;
        else
            nextIdxGlobal = elIdx(j)+1;
        end
        nextIdx = IEN(nextIdxGlobal);
        for i = 1:hDOF
              if hVecMesh(hDOF*(nextIdx-1)+i)>0 || gLogic_global(gDOF*nextIdxGlobal)
                  hIdx = hDOF*(elIdx(j)-1)+i;
                  hVec_global(hIdx) = mean([hVecCur(i),hVecMesh(hDOF*(nextIdx-1)+i)]);
              end
        end
    end
end
%output into struct
meshStruct.P = P;
meshStruct.IEN = IEN;
meshStruct.x = x;
meshStruct.y = y;
meshStruct.gVec = gVec_global;
meshStruct.hVec = hVec_global;
meshStruct.gNum = gNum;
meshStruct.axisymm = true;
meshStruct.numEl = numEl;
meshStruct.xVec = xVec;
meshStruct.yVec = yVec;
meshStruct.boundIndices = boundIndices;
meshStruct.nodesPerEl = 3;
end