function [rNew,zNew] = respaceNodes(r,z,transHeight,spacingVec)
% respace boundary nodes for a cell spreading on a flat surface. Boundary
% nodes for the bulk of the cell body are spaced by spacingVec(1) in terms 
% of arc length and near the surface (when z < transHeight) spacing is spacingVec(2). 
% There is a transition region in between where the spacing gradually increases from
% spacingVec(1) to spacingVec(2) (see code for details)
% INPUTS:
%   ** NOTE: boundary nodes should start at the top of the cell and go
%   clockwise around the cell, ending at the point (0,0) which is on the
%   substrate at the line of symmetry
%   r: r values for cell boundary nodes.
%   z: z values for cell boundary nodes
%   transHeight: z value (height off surface) where refined mesh starts
%   spacingVec: vector with two elements - first is the spacing in the main
%       cell body, second is spacing of refined mesh near surface
% OUTPUTS:
%   rNew: r values for newly spaced nodes
%   zNew: z values for newly spaced nodes
    
    % compute arc length (sParam)
    sParam = zeros(1,length(r));
    sParam(1) = 0;
    for i = 2:length(r)
        sParam(i) = sParam(i-1) + sqrt((r(i)-r(i-1)).^2 + (z(i)-z(i-1)).^2);
    end
    % find and eliminate possible repeat points
    repeatLogic = [diff(sParam) == 0,false];
    r = r(~repeatLogic);
    z = z(~repeatLogic);
    sParam = sParam(~repeatLogic);
    
    surfLogic = z<=eps; % eps is the smallest number for double precision
    startSurf = find(surfLogic,1,'first');
    if isempty(startSurf)
        startSurf = length(r); % only consider last node in contact in this case
    end
    
    startSpacing = spacingVec(1);
    minSpacing = spacingVec(2);
    if transHeight == 0 % then no extended refined boundary near surface
        sEndTrans = sParam(startSurf);
        rDense = [];
        zDense = [];
    else
        transIdx = find(z > transHeight, 1, 'last'); % last index before refined mesh
        if isempty(transIdx) || sParam(transIdx) < sParam(startSurf)/2 % then start refined mesh halfway
            [~,transIdx] = min(abs(sParam-sParam(startSurf)/2));
            sEndTrans = sParam(transIdx);
        else
            sEndTrans = interp1(z(transIdx:transIdx+1),sParam(transIdx:transIdx+1),transHeight);
        end
        denseSpan = sParam(startSurf)-sEndTrans;
        minSpacing = denseSpan/round(denseSpan/minSpacing);
        sDense = sEndTrans:minSpacing:sParam(startSurf);
        rDense = interp1(sParam,r,sDense);
        zDense = interp1(sParam,z,sDense);
    end
    % now compute position of "transition nodes", where spacing goes from
    % startSpacing to minSpacing. The spacing decreases by an increment of
    % ~minSpacing for each step (see documentation)
    m = round(startSpacing/minSpacing) + 1; % number of transition nodes
    transSpan = m*(m-1)*minSpacing / 2; % the distance spanned by this transition region 
    if transSpan > sEndTrans % then not enough room for full transition (shouldn't usually happen)
        transSpan = sEndTrans;
    end
    sTrans = sEndTrans - transSpan;
    % spacing gradually decreases as you approach the surface (see
    % documentation)
    for n = 1:m-2
        sTrans = [sTrans,sTrans(end)+(m-n)*minSpacing];  %#ok<AGROW>
    end
    if sTrans(1) < minSpacing % then there is no additional start region
        sTrans = sTrans(sTrans >= minSpacing);
        sTrans = [0,sTrans];
        sStartNew = [];
    else
        sStartTrans = sTrans(1);
        roundSpacing = sStartTrans/round(sStartTrans/startSpacing);
        sStartNew = 0:roundSpacing:sStartTrans-roundSpacing;
    end
    % now interpolate r and z values for new arc length values
    rTrans = interp1(sParam,r,sTrans);
    zTrans = interp1(sParam,z,sTrans);
    rStart = interp1(sParam,r,sStartNew);
    zStart = interp1(sParam,z,sStartNew);
    % now compute nodes on the substrate with spacing about equal to that
    % in the main cell body
    surfLength = sParam(end)-sParam(startSurf);
    surfSpacing = surfLength/round(surfLength/(startSpacing));
    sSurf = sParam(startSurf):surfSpacing:sParam(end);
    zSurf = interp1(sParam,z,sSurf);
    rSurf = interp1(sParam,r,sSurf);
    % concatenate all segments
    rNew = [rStart, rTrans, rDense, rSurf];
    zNew = [zStart, zTrans, zDense, zSurf];
end