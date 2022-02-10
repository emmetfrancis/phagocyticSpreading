function [rTest,zTest,elDissip,totDissip,splitDissip] = dissipationCalc(r,z,v_r,v_z)

    numCols = sum(z<=eps);
    numRows = length(r)/numCols;
    numEl = (numCols-1)*(numRows-1);
    rTest = zeros(1,numEl);
    zTest = zeros(1,numEl);
    elDissip = zeros(1,numEl);
    dissip1 = zeros(1,numEl);
    dissip2 = zeros(1,numEl);
    dissip3 = zeros(1,numEl);
    dissip4 = zeros(1,numEl);
    dissip5 = zeros(1,numEl);
    dissip6 = zeros(1,numEl);
    elVol = zeros(1,numEl);
    for i = 1:numEl
        curRow = floor((i-1) / (numCols - 1)) + 1;
        curCol = mod(i-1, (numCols - 1)) + 1;
        startingNode = numCols * (curRow - 1) + curCol;
        endingNode = numCols * curRow + curCol;
        globalIdx = [startingNode,startingNode+1,endingNode+1,endingNode];
        rCur = r(globalIdx);
        zCur = z(globalIdx);
        v_rCur = v_r(globalIdx);
        v_zCur = v_z(globalIdx);
        gaussPts = [-sqrt(3/5),0,sqrt(3/5)];
        gaussWts = [5/9,8/9,5/9];
        for n = 1:length(gaussPts)
            for m = 1:length(gaussPts)
                curEps = gaussPts(n);
                curEta = gaussPts(m);
                N_eps = [ -(1 + curEta) / 4, (1 + curEta) / 4, (1 - curEta) / 4, -(1 - curEta) / 4 ]';
                N_eta = [ (1 - curEps) / 4, (1 + curEps) / 4, -(1 + curEps) / 4, -(1 - curEps) / 4 ]';
                r_eps = sum(N_eps .* rCur);
                r_eta = sum(N_eta .* rCur);
                z_eps = sum(N_eps .* zCur);
                z_eta = sum(N_eta .* zCur);
                jVal = r_eps * z_eta - r_eta * z_eps;
                jMat = (1/jVal) * [z_eta, -r_eta; -z_eps, r_eps];
                NMat = horzcat(N_eps, N_eta);
                NMat = NMat * jMat; %  4x2, each row a is [dN_a/dx, dN_a/dy]
                N_cur = [ 0.25 * (1 - curEps) * (1 + curEta), 0.25 * (1 + curEps) * (1 + curEta),...
                    0.25 * (1 + curEps) * (1 - curEta), 0.25 * (1 - curEps) * (1 - curEta) ]';
                rVal = sum(N_cur .* rCur);
                zVal = sum(N_cur .* zCur);
                gWt = gaussWts(n)*gaussWts(m);
                rTest(i) = rTest(i) + gWt*2*pi*jVal*(rVal.^2);
                zTest(i) = zTest(i) + gWt*2*pi*jVal*rVal.*zVal;
                elVol(i) = elVol(i) + gWt*2*pi*jVal*rVal;
                dissip1(i) = dissip1(i) + 2*gWt*2*pi*jVal*rVal*sum(NMat(:,1).*v_rCur).^2;
                dissip2(i) = dissip2(i) + 2*gWt*2*pi*jVal*rVal*(sum(N_cur.*v_rCur)./rVal).^2;
                dissip3(i) = dissip3(i) + 2*gWt*2*pi*jVal*rVal*sum(NMat(:,2).*v_zCur).^2;
                dissip4(i) = dissip4(i) + gWt*2*pi*jVal*rVal*sum(NMat(:,2).*v_rCur).^2;
                dissip5(i) = dissip5(i) + 2*gWt*2*pi*jVal*rVal*sum(NMat(:,2).*v_rCur)*sum(NMat(:,1).*v_zCur);
                dissip6(i) = dissip6(i) + gWt*2*pi*jVal*rVal*sum(NMat(:,1).*v_zCur).^2;
%                 curDissip = (sum(NMat(:,2).*v_rCur) + sum(NMat(:,1).*v_zCur)).^2;
                elDissip(i) = elDissip(i) + dissip1(i) + dissip2(i) + dissip3(i) +...
                    dissip4(i) + dissip5(i) + dissip6(i);
%                 elDissip(i) = elDissip(i) + dissip5(i);
%                 elDissip(i) = elDissip(i) + gWt*2*pi*jVal*curDissip;
            end
        end
        rTest(i) = rTest(i)/elVol(i);
        zTest(i) = zTest(i)/elVol(i);
    end
    totDissip = sum(elDissip);
    elDissip = elDissip./elVol;
    splitDissip{1} = dissip1./elVol;
    splitDissip{2} = dissip2./elVol;
    splitDissip{3} = dissip3./elVol;
    splitDissip{4} = dissip4./elVol;
    splitDissip{5} = dissip5./elVol;
    splitDissip{6} = dissip6./elVol;
%     vzDissip = sum(vzDissip);
%     vrDissip = sum(vrDissip);
end