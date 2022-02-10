%% Load video data
loadMore = true;
n = 0;
videoData = struct;
while loadMore
    n = n+1;
    curFolder = uigetdir('',sprintf('Choose folder number %d',n));
    if curFolder == 0
        loadMore = false;
        break
    end
    videoData(n).folder = curFolder;
    try loadStruct = load(fullfile(curFolder,'dataStruct.mat'));
        dataStruct = loadStruct.dataStruct;
        videoData(n).rStored = dataStruct.rStored;
        videoData(n).zStored = dataStruct.zStored;
        videoData(n).v_rStored = dataStruct.v_rStored;
        videoData(n).v_zStored = dataStruct.v_zStored;
        videoData(n).pStored = dataStruct.pStored;
        videoData(n).timeVec = dataStruct.timeVec;
        videoData(n).IENStored = dataStruct.IENStored;
        videoData(n).boundStored = dataStruct.boundStored;
        videoData(n).theta_nStored = dataStruct.theta_nStored;
    catch
        rData = load(fullfile(curFolder,'rData.mat'));
        zData = load(fullfile(curFolder,'zData.mat'));
        v_rData = load(fullfile(curFolder,'v_rData.mat'));
        v_zData = load(fullfile(curFolder,'v_zData.mat'));
        pData = load(fullfile(curFolder,'pData.mat'));
        timeData = load(fullfile(curFolder,'timeVec.mat'));
        try theta_nData = load(fullfile(curFolder,'theta_nData.mat'));
            videoData(n).theta_nStored = theta_nData.theta_nStored;
        catch % single phase
            videoData(n).theta_nStored = rData.rStored;
            for i = 1:length(rData.rStored)
                videoData(n).theta_nStored{i} = .001 * ones(size(rData.rStored{i}));
            end
        end
        videoData(n).rStored = rData.rStored;
        videoData(n).zStored = zData.zStored;
        videoData(n).v_rStored = v_rData.v_rStored;
        videoData(n).v_zStored = v_zData.v_zStored;
        videoData(n).pStored = pData.pStored;
        videoData(n).timeVec = timeData.timeVec;
    end
end
NDlg = inputdlg('Enter N values in order with commas');
NTxt = NDlg{1};
NVals = textscan(NTxt,'%.6f,');
NVals = NVals{1};
% maxTimes = [200,500,1000,1000,1000];
% maxTimes = [200,200,200,200];
% startTimes = [42,42,42,42];
maxTimes = 200*ones(1,length(NVals));
startTimes = zeros(1,length(NVals)) + 42;

%% View insideFlow videos
figure
skipVal = [5,10];
saveVid = true;
discreteLogic = true;
CaLogic = false;
CaStrings = {'Standard Conditions','2 \muM BAPTA-AM'};
if saveVid
    frameRate = 30;
    path = uigetdir('Choose where to save video');
    myVideo = VideoWriter(fullfile(path,'Movie5.mp4'),'MPEG-4');
    myVideo.FrameRate = frameRate;
    open(myVideo)
end
set(gcf, 'Position', [500 200 800 500])

for k = 1:length(videoData)
    rStored = videoData(k).rStored;
    zStored = videoData(k).zStored;
    v_rStored = videoData(k).v_rStored;
    v_zStored = videoData(k).v_zStored;
    pStored = videoData(k).pStored;
    timeVec = videoData(k).timeVec;
    theta_nStored = videoData(k).theta_nStored;
%     if k <= 2
%         timeVec = timeVec * 1660/600;
%     end
    timeSample = startTimes(k):.5:maxTimes(k);
%     timeSample = 0:.2:20;
    idxSample = zeros(1,length(timeSample));
    for i = 1:length(timeSample)
        [~,idxSample(i)] = min(abs(timeVec-timeSample(i)));
    end
    vNorm = 0.5;
    terms = 'all';
    framesCompleted = false(1,length(timeSample));
    % initial mesh is cartesian with edge points
    startIdx = 1;
    rStart = rStored{startIdx};
    zStart = zStored{startIdx};
    spacing = .033333;
%     spacing = .0333333/2;
    rTest = -2:spacing:2;
    zTest = 0:2.5*spacing:1;
    [rMesh,zMesh] = meshgrid(rTest,zTest);
    for i = 2:2:size(rMesh,1)
        rMesh(i,:) = rMesh(i,:)+spacing/2;
    end
    rAll = rMesh(:)';
    zAll = zMesh(:)';
    if length(rStart) < 200
        rEdge = [rStart,-rStart(end-1:-1:2)];
        zEdge = [zStart,zStart(end-1:-1:2)];
    else
        if isfield(videoData(k),'boundStored')
            boundIndices = videoData(k).boundStored{startIdx};
        else
            numCols = sum(zStart<eps);
            numRows = length(zStart)/numCols;
            topIndices = 1:numCols-1;
            rightIndices = numCols*(1:numRows-1);
            bottomIndices = numRows*numCols:-1:((numRows-1)*numCols+2);
            leftIndices = numCols*(numRows-1:-1:1)+1;
            boundIndices = [topIndices,rightIndices,bottomIndices,leftIndices(1)];
        end
        rEdge = [rStart(boundIndices);-rStart(boundIndices(end-1:-1:1))];
        zEdge = [zStart(boundIndices); zStart(boundIndices(end-1:-1:1))];
    end
    [in,on] = inpolygon(rAll,zAll,rEdge,zEdge);
    rAll = rAll(in);
    zAll = zAll(in);
    v_zPrev = 0;
    v_rPrev = 0;
    rMeshOrig = rMesh;
    zMeshOrig = zMesh;
    
    
    
    
    for i = startIdx:length(rStored)
        if all(framesCompleted) || timeVec(i) > timeSample(end) + 5 || i > length(v_rStored)
            break
        end
        curTime = timeVec(i);
        rCur = rStored{i};
        zCur = zStored{i};
        if isempty(rCur)
            break
        end
        
        cla
        x = [-15, -15, 15, 15];
        y = [0, -.5, -.5, 0];
        patch('XData',x,'YData',y,'FaceColor',[.8 .8 .8])
        hold on
        if length(rCur) < 200
            curCoeff = coeffStored{i+1};
            patch('XData',8.5*[rCur,-rCur(end-1:-1:2)],'YData',8.5*[zCur,zCur(end-1:-1:2)],...
                'FaceColor',[.5 .5 .5],'LineWidth',2,'EdgeColor','k')
            sParam = zeros(1,length(rCur));
            for j = 2:length(rCur)
                sParam(j) = sParam(j-1) + sqrt((rCur(j)-rCur(j-1))^2 + (zCur(j)-zCur(j-1))^2);
            end
            [~,fineIdx] = min(diff(diff(sParam)));
            fineIdx = fineIdx + 2;
            contactIdx = find(zCur <= 0, 1, 'first');
            skipVal = round(mean(diff(sParam(1:fineIdx-1))) / mean(diff(sParam(fineIdx:contactIdx))));
            edgeIdx = [1:3:fineIdx,fineIdx + skipVal:3*skipVal:contactIdx-1];
            rEdge = [rCur(edgeIdx), -rCur(edgeIdx(end:-1:2))];
            zEdge = [zCur(edgeIdx), zCur(edgeIdx(end:-1:2))];
            % first plot edge velocities
            [v_rEdge,v_zEdge,~] = meshFlowCalc(rEdge,zEdge,curCoeff,terms);
            v_rEdge(rEdge < 0) = -v_rEdge(rEdge < 0);
            %     quiver(rEdge*8.5,zEdge*8.5,v_rEdge/vNorm,v_zEdge/vNorm,'LineWidth',1,'Color','k','AutoScale','off','MaxHeadSize',.1)
            % now plot inner velocities
            [v_r,v_z,~] = meshFlowCalc(rAll,zAll,curCoeff,terms);
            v_r(rAll < 0) = -v_r(rAll < 0);
            %         v_z(zAll < .02) = 0;
            %         v_r(zAll < .02) = 0;
            quiver(rAll*8.5,zAll*8.5,v_r/vNorm,v_z/vNorm,'LineWidth',1,...
                'Color','b','AutoScale','off','MaxHeadSize',.1)
        else
            numCols = sum(zCur<eps);
            numRows = length(zCur)/numCols;
            v_rStored{i} = v_rStored{i}*.01/600;
            v_zStored{i} = v_zStored{i}*.01/600;
            v_rCur = v_rStored{i};
            v_zCur = v_zStored{i};
            pCur = pStored{i}-100;
            theta_nCur = theta_nStored{i};
            if isfield(videoData(k),'boundStored')
                boundIndices = videoData(k).boundStored{i};
            else
                topIndices = 1:numCols-1;
                rightIndices = numCols*(1:numRows-1);
                bottomIndices = numRows*numCols:-1:((numRows-1)*numCols+2);
                leftIndices = numCols*(numRows-1:-1:1)+1;
                boundIndices = [topIndices,rightIndices,bottomIndices,leftIndices(1)];
            end
            
            rEdge = [rCur(boundIndices);-rCur(boundIndices(end-1:-1:1))];
            zEdge = [zCur(boundIndices); zCur(boundIndices(end-1:-1:1))];
            if ~CaLogic
                patch('XData',8.5*rEdge,'YData',8.5*zEdge,...
                    'FaceColor',[.85 .85 .85],'LineWidth',2,'EdgeColor','k')
            else
                MFIVal = interp1(CaData{k}(:,1)+20,CaData{k}(:,2),curTime);
                patch('XData',8.5*rEdge,'YData',8.5*zEdge,...
                'FaceColor',[0 .1+.9*MFIVal/5.5 0],'LineWidth',2,'EdgeColor','k')
                cmap = [zeros(256,1),.1+.9*(0:255)'/255,zeros(256,1)];
                colormap(cmap)
                colorbar(gca,'Ticks',[1/5.5,3/5.5,5/5.5],'TickLabels',...
                    {'[Ca^{2+}]_0','3\times[Ca^{2+}]_0','5\times[Ca^{2+}]_0'})
            end
%             patch('XData',8.5*rEdge,'YData',8.5*zEdge,...
%                  'FaceColor','none','LineWidth',2,'EdgeColor','k')
%             rMeshCur = rMeshOrig(:)';
%             zMeshCur = zMeshOrig(:)';
%             [in,on] = inpolygon(rMeshCur,zMeshCur,rEdge,zEdge);
%             rMeshCur = rMeshCur(in | on);
%             zMeshCur = zMeshCur(in | on);
            v_rInterp = scatteredInterpolant([rCur;-rCur],[zCur;zCur],[v_rStored{i};-v_rStored{i}]);
            v_r = v_rInterp(rAll,zAll);
%             v_r = v_rInterp(rMeshCur,zMeshCur);
            v_rEdge = v_rInterp(rEdge,zEdge);
            v_zInterp = scatteredInterpolant([rCur;-rCur],[zCur;zCur],[v_zStored{i};v_zStored{i}]);
            v_z = v_zInterp(rAll,zAll);
%             v_z = v_zInterp(rMeshCur,zMeshCur);
            v_zEdge = v_zInterp(rEdge,zEdge);
%             if i > startIdx
%                 v_rPrev = v_rInterpPrev(rAll,zAll);
%                 v_zPrev = v_zInterpPrev(rAll,zAll);
%                 v_rChange = abs((v_r-v_rPrev)./v_rPrev);
%                 v_zChange = abs((v_z-v_zPrev)./v_zPrev);
%                 v_rOL = v_rChange > 10;
%                 v_zOL = v_zChange > 10;
% %                 if any(v_rOL)
% %                     v_r(v_rOL) = v_rPrev(v_rOL);
% %                     v_rInterp = scatteredInterpolant(rAll',zAll',v_r');
% %                 end
%                 if any(v_zOL)
%                     v_z(v_zOL) = v_zPrev(v_zOL);
%                     v_zInterp = scatteredInterpolant(rAll',zAll',v_z');
%                 end
%         end
%             v_rMesh = reshape(v_rCur,[numCols,numRows])';
%             v_zMesh = reshape(v_zCur,[numCols,numRows])';
%             [vr_r,vr_z] = gradient(v_rMesh);
%             vrGradMag = sqrt(vr_r.^2 + vr_z.^2) / mean(abs(v_rCur));
%             [vz_r,vz_z] = gradient(v_zMesh);
%             vzGradMag = sqrt(vz_r.^2 + vz_z.^2) / mean(abs(v_zCur));
            
%             topLogic = zAll > 0.9*max(zAll);
% %             v_zSpike = vzGradMag(:) > 2;
% %             v_rSpike = vrGradMag(:) > 10;
%             v_zSpike = topLogic & v_z > -0.1*min(v_z);
%             if any(v_zSpike)
%                 rAll = rAll + v_r*(.025/2023) * (timeVec(i+1)-timeVec(i))/8.5e-4;
%                 zAll = zAll + v_z*(.025/2023) * (timeVec(i+1)-timeVec(i))/8.5e-4;
%                 zAll(zAll < 0) = 0;
%                 [in,on] = inpolygon(rAll,zAll,rEdge,zEdge);
%                 rAll = rAll(in);
%                 zAll = zAll(in);
%                 continue
%             end
            vIndNorm = sqrt(v_r.^2 + v_z.^2);
%             quiver(rAll*8.5,zAll*8.5,v_r/vNorm,v_z/vNorm,'LineWidth',1,...
%                 'Color','b','AutoScale','off','MaxHeadSize',.1)
%             quiver(rMeshCur*8.5,zMeshCur*8.5,v_r/vNorm,v_z./vNorm,'LineWidth',1,...
%                 'Color','b','AutoScale','off','MaxHeadSize',.1)
%             plot(rAll*8.5,zAll*8.5,'Marker','o','MarkerSize',2,'LineStyle','none',...
%                 'MarkerEdgeColor','b','MarkerFaceColor','b')
            
%             pMesh = reshape(pCur,[numCols,numRows])';
%             surf(rMesh*8.5,zMesh*8.5,pMesh,pMesh,'FaceColor','interp','EdgeColor','none')
%             hold on
%             surf(-rMesh*8.5,zMesh*8.5,pMesh,pMesh,'FaceColor','interp','EdgeColor','none')
%             caxis([-120 0])
%             view([0 90])
%             eIdx = [1:3:length(boundIndices),length(rEdge):-3:length(boundIndices)+1];
%             quiver(rEdge(eIdx)*8.5,zEdge(eIdx)*8.5,v_rEdge(eIdx)/vNorm,v_zEdge(eIdx)/vNorm,'LineWidth',1,...
%                 'Color','k','AutoScale','off','MaxHeadSize',.1)
%             theta_nCur(theta_nCur < 0.001) = 0.001;
%             if isfield(videoData(k),'IENStored')
%                 IEN = videoData(k).IENStored{i};
%                 IEN = reshape(IEN,[3,length(IEN)/3])';
%                 trisurf(IEN,8.5*rCur,8.5*zCur,theta_nCur,theta_nCur,'FaceColor','interp','EdgeColor','none')
%                 trisurf(IEN,-8.5*rCur,8.5*zCur,theta_nCur,theta_nCur,'FaceColor','interp','EdgeColor','none')
%             else
%                 rMesh = reshape(rCur,[numCols,numRows])';
%                 zMesh = reshape(zCur,[numCols,numRows])';
%                 theta_nMesh = reshape(theta_nCur,[numCols,numRows])';
%                 surf(rMesh*8.5,zMesh*8.5,theta_nMesh,theta_nMesh,'FaceColor','interp','EdgeColor','none')
%                 hold on
%                 surf(-rMesh*8.5,zMesh*8.5,theta_nMesh,theta_nMesh,'FaceColor','interp','EdgeColor','none')
%             end
%             caxis([.00 .015])
%             view([0 90])
%             h = colorbar;
%             set(get(h,'title'),'string','\theta_n');
%             plot(8.5*rEdge,8.5*zEdge,'-k','LineWidth',2)
        end
        if discreteLogic
            NCur = NVals(k);
%             NCur = 30;
            NSpacing = sqrt(1/NCur);
            %     NShift = mod(0,sqrt(1/NCur));
            fractionVal = 0.0;
            NShift = fractionVal * NSpacing;
            NLoc = NShift:NSpacing:20;
            NAll = [-NLoc(end:-1:1),NLoc];
            curRad = rEdge(find(zEdge<=eps,1,'first'))*8.5;
            if NCur <= 1000
                plot(NAll(abs(NAll) <= curRad+1e-3),zeros(1,sum(abs(NAll)<=curRad+1e-3)),'ok','MarkerFaceColor','k','MarkerSize',5)
                plot(NAll(abs(NAll) > curRad+1e-3),zeros(1,sum(abs(NAll)>curRad+1e-3)),'ok','MarkerFaceColor','w','MarkerSize',5)
            else
                NBoundIdx = find(NAll <= curRad, 1, 'last');
                patch('XData',[0,NAll(NBoundIdx),NAll(NBoundIdx),0],...
                    'YData',[.02,.02,-.02,-.02],'FaceColor','k','EdgeColor','k')
                patch('XData',[NAll(NBoundIdx),NAll(end),NAll(end),NAll(NBoundIdx)],...
                    'YData',[.02,.02,-.02,-.02],'FaceColor','w','EdgeColor','k')
            end
        end
        daspect([1 1 1])
%         xlim([4,10.5])
%         ylim([-.1,3])
        xlim([5,10])
        ylim([-.1,3])
        xlabel('r [\mum]')
        ylabel('z [\mum]')
        grid off
        prettyGraph
        drawnow
        hold off
        % update rAll and zAll
        rAll = rAll + v_r*(.025/2023) * (timeVec(i+1)-timeVec(i))/8.5e-4;
        zAll = zAll + v_z*(.025/2023) * (timeVec(i+1)-timeVec(i))/8.5e-4;
        rAll = rAll + v_r * (timeVec(i+1)-timeVec(i))/8.5;
        zAll = zAll + v_z * (timeVec(i+1)-timeVec(i))/8.5;
        zAll(zAll < 0) = 0;
        [in,on] = inpolygon(rAll,zAll,rEdge,zEdge);
        rAll = rAll(in);
        zAll = zAll(in);
        
        curIdx = find(idxSample == i);
        if length(curIdx) > 1
            curIdx = curIdx(1);
        end
        if ~isempty(curIdx)
%             title(sprintf('Adhesion-driven spreading, \\rho_{l} = %.1d %%', NVals(k)/100), ...
%                 sprintf('t = %.2f s',timeSample(curIdx)))
%             title(sprintf('Protrusion-only model'), ...
%                 sprintf('t = %.2f s',timeSample(curIdx)))
%               title(sprintf('Protrusion-only model, s_0 = %.1f \\mum', NVals(k)), ...
%                 sprintf('t = %.2f s',timeSample(curIdx)))
%             title(sprintf('Protrusion-driven spreading, \\rho_{l} = %.1d %%', NVals(k)/100), ...
%                 sprintf('t = %.2f s',timeSample(curIdx)))
            title(sprintf('Discrete model, \\rho_{l} = %d \\mum^{-2}', NVals(k)), ...
                sprintf('t = %.2f s',timeSample(curIdx)))
%               title(sprintf('%s', CaStrings{k}), ...
%                 sprintf('t = %.2f s',timeSample(curIdx)))
%             title('Network volume fraction', ...
%                 sprintf('t = %.2f s',timeSample(curIdx)))
            if saveVid
%                 img = getframe(gcf);
%                 img = img.cdata;
                img = print('-RGBImage');
%                 img = print('-RGBImage','-r80');
                writeVideo(myVideo, img)
            end
            framesCompleted(curIdx) = true;
        end
    end
    if saveVid
        for i = 1:frameRate
%             img = getframe(gcf);
%             img = img.cdata;
            img = print('-RGBImage');
%             img = print('-RGBImage','-r80');
            writeVideo(myVideo, img)
        end
    end
    
end

if saveVid
    close(myVideo)
end