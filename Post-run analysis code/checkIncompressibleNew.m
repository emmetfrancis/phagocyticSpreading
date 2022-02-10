function [rTest,zTest,elDil,totDil,mse] = checkIncompressibleNew(r,z,v_r,v_z)

    numCols = sum(z<=eps);
    numRows = length(r)/numCols;
    numEl = (numCols-1)*(numRows-1);
    rTest = zeros(1,numEl);
    zTest = zeros(1,numEl);
    elDil = zeros(1,numEl);
    elDilSq = zeros(1,numEl);
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
                NMat(:,1) = NMat(:,1) + N_cur/rVal;
                gWt = gaussWts(n)*gaussWts(m);
                elVol(i) = elVol(i) + gWt*2*pi*jVal*rVal;
                vDiv = sum(NMat(:,1).*v_rCur) + sum(NMat(:,2).*v_zCur);
                elDil(i) = elDil(i) + gWt*2*pi*jVal*rVal*vDiv;
                elDilSq(i) = elDilSq(i) + gWt*2*pi*jVal*rVal*vDiv^2;
            end
        end
        % record centroid to id element
        curPolygon = polyshape(rCur,zCur);
        [rC,zC] = centroid(curPolygon);
        rTest(i) = rC;
        zTest(i) = zC;
    end
    totDil = sum(elDil) / sum(elVol);
    mse = sum(elDilSq) / sum(elVol);
    elDil = elDil ./ elVol;
end