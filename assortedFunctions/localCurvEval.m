function [mCurv,pCurv,phiVals] = localCurvEval(r,z,inclDiscont)

    fitType = 3; % cubic fit
    numFitIndices = 9;
    halfInt = floor(numFitIndices/2);
    sParam = zeros(1,length(r));
    for i = 2:length(sParam)
        sParam(i) = sParam(i-1) + sqrt((r(i)-r(i-1))^2 + (z(i)-z(i-1))^2);
    end
    surfLogic = z <= eps;
    startSurf = find(surfLogic, 1, 'first');
    
    r360 = horzcat(-r(halfInt+1:-1:2), r); % add in symmetric points
    z360 = horzcat(z(halfInt+1:-1:2), z);
    surfLogic360 = z360 <= eps;
    startSurf360 = find(surfLogic360, 1, 'first');
    
    [mCurv,pCurv,phiVals] = localCurv(r360(1:startSurf360),...
        z360(1:startSurf360),fitType,numFitIndices);
    phiValClose = phiVals(end-halfInt);
    mCurvClose = mCurv(end-halfInt);
    pCurvFix = isinf(pCurv) | isnan(pCurv);
    pCurv(pCurvFix) = mCurv(pCurvFix);
    
    % set constant curvature near surface and zero on surface
    if inclDiscont
        phiIncrement = (pi-phiValClose) / halfInt;
        sIncrement = mean(diff(sParam(startSurf-halfInt:startSurf)));
        phiValsSurf = [phiValClose:phiIncrement:pi-phiIncrement, pi*ones(1,sum(surfLogic))];
        mCurvSurf = [(phiIncrement/sIncrement)*ones(1,halfInt), zeros(1,sum(surfLogic))];
        pCurvSurf = [sin(phiValsSurf(1:halfInt))./r(startSurf-halfInt:startSurf-1), zeros(1,sum(surfLogic))];
    else
        linFit_r = polyfit(sParam(startSurf-3:startSurf),r(startSurf-3:startSurf),1);
        linFit_z = polyfit(sParam(startSurf-3:startSurf),z(startSurf-3:startSurf),1);
        phiContact = atan2(-linFit_z(1),linFit_r(1));
        if mCurvClose < 1e-6
            phiValsSurf  = [phiContact*ones(1,halfInt), pi*ones(1,sum(surfLogic))];
            mCurvSurf = [zeros(1,halfInt), zeros(1,sum(surfLogic))];
            pCurvSurf = [zeros(1,halfInt), zeros(1,sum(surfLogic))];
        else
            %         phiSteps = mCurvClose * diff(sParam(startSurf-halfInt:startSurf));
            %         phiApproach = phiValClose * ones(1,halfInt);
            %         for i = 2:halfInt
            %             phiApproach(i) = phiApproach(i-1) + phiSteps(i-1);
            %         end
            rCenter = r(startSurf) - (1/mCurvClose)*sin(phiContact);
            zCenter = -(1/mCurvClose)*cos(phiContact);
            phiIncr = mCurvClose * (sParam(startSurf)-sParam(startSurf-1)); % continuing at smallest spacing
            sIncrStart = sParam(startSurf-halfInt+1) - sParam(startSurf-halfInt);
            numImagPts = ceil(halfInt*sIncrStart / (phiIncr/mCurvClose));
            phiImag = phiContact + (1:numImagPts) * phiIncr;
            rImag = (1/mCurvClose)*sin(phiImag) + rCenter;
            zImag = (1/mCurvClose)*cos(phiImag) + zCenter;
            [mCurvApproach,pCurvApproach,phiApproach] = localCurv([r(startSurf-numFitIndices:startSurf),rImag],...
                [z(startSurf-numFitIndices:startSurf),zImag],fitType,numFitIndices);
            phiValsSurf  = [phiApproach(halfInt+1:2*halfInt), pi*ones(1,sum(surfLogic))];
            mCurvSurf = [mCurvApproach(halfInt+1:2*halfInt), zeros(1,sum(surfLogic))];
            pCurvSurf = [pCurvApproach(halfInt+1:2*halfInt), zeros(1,sum(surfLogic))];
        end
%         phiValsSurf = [phiApproach, pi*ones(1,sum(surfLogic))];
%         mCurvSurf = [mCurvClose*ones(1,halfInt), zeros(1,sum(surfLogic))];
%         pCurvSurf = [sin(phiValsSurf(1:halfInt))./r(startSurf-halfInt:startSurf-1), zeros(1,sum(surfLogic))];
    end
    
    mCurv = [mCurv(halfInt+1:end-halfInt-1),mCurvSurf];
    pCurv = [pCurv(halfInt+1:end-halfInt-1),pCurvSurf];
    phiVals = [phiVals(halfInt+1:end-halfInt-1),phiValsSurf];
    
end