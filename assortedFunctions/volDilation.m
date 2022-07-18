function [rNew,zNew] = volDilation(r,z,targetVol,volTol,topPlateLogic)
    [~,vol] = SAVolCalc(r,z,'trapz');
    fitType = 3; % cubic fit
    rOrig = r;
    zOrig = z;
    numFitIndices = 11; % fit cur point and two to left and right
    halfInt = floor(numFitIndices/2);
    r360 = horzcat(-r(numFitIndices+1:-1:2), r, -r(end-1:-1:end-numFitIndices)); % add in symmetric points
    z360 = horzcat(z(numFitIndices+1:-1:2), z, z(end-1:-1:end-numFitIndices));
    sParam = zeros(1,length(r));
    for i = 2:length(sParam)
        sParam(i) = sParam(i-1) + sqrt((r(i)-r(i-1))^2 + (z(i)-z(i-1))^2);
    end
    startingGap = mean(diff(sParam(1:10)));
    startDense360 = find(abs(diff(sParam)-startingGap) > .001,1,'first');
    surfLogic360 = z360 <= 1e-3;
    startSurf360 = find(surfLogic360, 1, 'first');
    surfLogic = z <= 1e-3;
    startSurf = find(surfLogic, 1, 'first');
    int1 = 1:startDense360+numFitIndices;
    int2 = startDense360-numFitIndices:startSurf360;
    if isempty(startDense360) || (startDense360-halfInt-1) < numFitIndices+1 % all points set to same density
        phiVals2 = [];
        int1 = 1:startSurf360;
        [~,~,phiVals1] = polyCurv([],r360(int1),z360(int1),fitType,numFitIndices);
        phiValsMid = [];
        phiValClose = phiVals1(end-halfInt);
    else
        [~,~,phiVals1] = polyCurv([],r360(int1),z360(int1),fitType,numFitIndices-4); %only need ~7 to fit this
        [~,~,phiVals2] = polyCurv([],r360(int2),z360(int2),fitType,numFitIndices);
        phiValsMid = phiVals2(halfInt+2:halfInt+numFitIndices+1);
        phiValClose = phiVals2(end-halfInt);
    end
%     [mCurvSurf,pCurvSurf,phiValsSurf] = polyCurv([],r360(surfInt),z360(surfInt),fitType,numFitIndices);

    % set constant curvature near surface and zero on surface
    phiIncrement = (pi-phiValClose) / halfInt;
    phiValsSurf = [phiValClose:phiIncrement:pi, pi*ones(1,length(r)-startSurf)];
    
    if isempty(startDense360) || (startDense360-halfInt-1) < numFitIndices+1
        phiVals = horzcat(phiVals1(numFitIndices+1:end-halfInt-1),phiValsSurf);
    else
        phiVals = horzcat(phiVals1(numFitIndices+1:startDense360-halfInt-1),phiValsMid,...
            phiVals2(halfInt+numFitIndices+2:end-halfInt-1),phiValsSurf);
    end
    
    % calculate amount dilation
    volCorrect = targetVol - vol;
    trialStep = 10*volCorrect*(startingGap/vol);
    trialStep = trialStep*ones(1,sum(~surfLogic));
    dummyVar = 0:20;
    trialStep(end-20:end) = trialStep(end-20:end) .* exp(-dummyVar/4);
    if topPlateLogic
        endTop = find(z == max(z),1,'last');
        trialStep(1:endTop) = 0;
        trialStep(endTop+21:-1:endTop+1) = trialStep(endTop+20:-1:endTop) .* exp(-dummyVar/4);
    end
    numIt = 0;
    while abs(vol-targetVol) > (volTol/100)*targetVol
        volOld = vol;
        rOld = r;
        zOld = z;
        r(~surfLogic) = r(~surfLogic) + trialStep.*sin(phiVals(~surfLogic));
        z(~surfLogic) = z(~surfLogic) + trialStep.*cos(phiVals(~surfLogic));
        r(r<0) = rOld(r<0);
        z(z<0) = zOld(z<0);
        [~,vol] = SAVolCalc(r,z,'trapz');
        closeness = (targetVol-vol) / (targetVol-volOld);
        trialStep = trialStep*closeness;
        numIt = numIt + 1;
        if numIt > 100
            rNew = rOrig;
            zNew = zOrig;
            return
        end
    end
    z(z<0) = 1e-12;
    rNew = r;
    zNew = z;
end
    