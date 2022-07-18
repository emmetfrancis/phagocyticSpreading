contactAngle = pi/3;
R0 = 8.5/2;
RCur = (4*R0^3 / ((1+cos(contactAngle))^2 * (2-cos(contactAngle))))^(1/3);
h = RCur * (1+cos(contactAngle));
sSurf = (pi-contactAngle)*RCur;
sEndTrans = sSurf;
sStartTrans = sSurf - 1;
transSpan = sEndTrans-sStartTrans;
startSpacing = .02*8.5;
minSpacing = .2*.05;
% m = round(2*transSpan / (minSpacing + startSpacing));
m = round(startSpacing/minSpacing) +  1;
transSpan = m*(m-1)*minSpacing/2;
% sInc = (startSpacing - minSpacing) / (m-1);
sTrans = minSpacing;
n = 0;
while sTrans(end) < transSpan-minSpacing
    n = n+1;
%     sTrans = [sTrans,sTrans(end)+n*sInc+minSpacing];
    sTrans = [sTrans,sTrans(end)+n*minSpacing];
end
sTrans = sEndTrans-sTrans;
sTrans = sTrans(end:-1:1);
sStartTrans = sTrans(1);
roundSpacing = sStartTrans/round(sStartTrans/startSpacing);
sStartNew = 0:roundSpacing:sStartTrans-roundSpacing;
angSample = (1/RCur)*[sStartNew,sTrans];
rBody = RCur*sin(angSample);
zBody = RCur*cos(angSample) + (h-RCur);
rContact = RCur*sin(contactAngle);
rSurf = rContact:-rContact/10:0;
zSurf = zeros(1,length(rSurf));
r = [rBody,rSurf];
z = [zBody,zSurf];
sParam = zeros(1,length(r));
for i = 2:length(sParam)
    sParam(i) = sParam(i-1) + sqrt((r(i)-r(i-1))^2 + (z(i)-z(i-1))^2);
end
%     [~,fineIdx] = min(diff(diff(sParam)));
%     fineIdx = fineIdx + 2;
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
% make mesh for velocities
gSides_v = [ 0,0; 0,0; 0,0; 0,0 ]; %top,right,bottom,left
gSidesLogic_v = [0,0; 0,0; 1,1; 0,0]; % left is axis of symm here
hSides_v = [ 1,1; 1,1; 0,0; 0,0 ];
hSidesLogic_v = [ 1,1; 1,1; 0,0; 1,1 ];
gCell_v = {gSides_v,gSidesLogic_v};
hCell_v = {hSides_v,hSidesLogic_v};
vMeshStruct = makeMeshNew(shapeCell,gCell_v,hCell_v);
xVec = vMeshStruct.xVec;
yVec = vMeshStruct.yVec;
figure
FEMPlot_v(xVec,yVec,zeros(size(xVec)),zeros(size(xVec)),numRows)
FEMPlot_v(-xVec,yVec,zeros(size(xVec)),zeros(size(xVec)),numRows)
xSurf = [-5,-5,5,5];
ySurf = [0, -.05, -.05, 0];
patch('XData',xSurf,'YData',ySurf,'FaceColor','k')
% xlim([3.5 4.5])
% ylim([-.05 1])
%%
figure
[mCurv,pCurv,phiVals] = localCurvEval(r,z,false);
plot(r,z,'LineWidth',2,'Color','k')
hold on
for i = 1:surfIdx-1
    plot([r(i),r(i)-.2*sin(phiVals(i))],[z(i),z(i)-.2*cos(phiVals(i))],'k','LineWidth',1)
end
daspect([1 1 1])
% xlim([3.5 4.5])
% ylim([0 .5])
prettyGraph

%% test alt spacing
startSpacing = .02*8.5;
minSpacing = .05 * .05;
r0Adh = .05;
z0 = -.05;
f0 = 400;
N = 1000;
maxAdhDist = 10;
figure
m = round(startSpacing/minSpacing) +  1;
transSpan = m*(m-1)*minSpacing/2;
% sInc = (startSpacing - minSpacing) / (m-1);
sTrans = minSpacing;
n = 0;
if transSpan > sEndTrans
    transSpan = sEndTrans;
end
while sTrans(end) < transSpan-minSpacing
    n = n+1;
    %     sTrans = [sTrans,sTrans(end)+n*sInc+minSpacing];
    sTrans = [sTrans,sTrans(end)+n*minSpacing];
end
sTrans = sEndTrans-sTrans;
sTrans = sTrans(end:-1:1);
sStartTrans = sTrans(1);
roundSpacing = sStartTrans/round(sStartTrans/startSpacing);
sStartNew = 0:roundSpacing:sStartTrans-roundSpacing;
angSample = (1/RCur)*[sStartNew,sTrans];
rBody = RCur*sin(angSample);
zBody = RCur*cos(angSample) + (h-RCur);
rContact = RCur*sin(contactAngle);
rSurf = rContact:-rContact/10:0;
zSurf = zeros(1,length(rSurf));
r = [rBody,rSurf];
z = [zBody,zSurf];
sParam = zeros(1,length(r));
for i = 2:length(sParam)
    sParam(i) = sParam(i-1) + sqrt((r(i)-r(i-1))^2 + (z(i)-z(i-1))^2);
end
adhStress_r = zeros(1,length(r));
adhStress_z = zeros(1,length(z));
surfLogic = z <= eps;
rEdge = r(find(surfLogic,1,'first'));
adhForce = adhForcePointDoubleTrapz_r4(r(~surfLogic),z(~surfLogic)-z0,maxAdhDist*10,r0Adh,f0,N,rEdge,N);
% not considering adhesion forces once on surface, so set to zero
adhStressOnly_r = [adhForce(1,:), zeros(1,sum(surfLogic))];
adhStressOnly_z = [adhForce(2,:), zeros(1,sum(surfLogic))];
figure
% plot(sParam,sqrt(adhStressOnly_r.^2 + adhStressOnly_z.^2))
hold on
% test different refinements...
sInt = (.1:.1:1) * minSpacing;
for j = 1:length(sInt)
    for i = 2:size(adhForce,2)
        sIntVals = sParam(i)-sInt(j):sInt(j):sParam(i)+sInt(j);
        rIntVals = interp1(sParam,r,sIntVals);
        rStressVals = interp1(sParam,adhStressOnly_r,sIntVals);
        zStressVals = interp1(sParam,adhStressOnly_z,sIntVals);
        rAdhInt = trapz(sIntVals,rStressVals .* rIntVals);
        zAdhInt = trapz(sIntVals,zStressVals .* rIntVals);
        areaInt = trapz(sIntVals,rIntVals);
        adhStress_r(i) = rAdhInt / areaInt;
        adhStress_z(i) = zAdhInt / areaInt;
    end
    curError = sqrt((adhStress_r-adhStressOnly_r).^2 + (adhStress_z-adhStressOnly_z).^2);
    plot(sParam,curError)
%     plot(sParam,sqrt(adhStress_r.^2 + adhStress_z.^2))
    hold on
end