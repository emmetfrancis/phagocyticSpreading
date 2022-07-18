function [gCell,hCell] = setCellBC(r,z,condition)

surfLogic = z <= eps;
startSurf = find(surfLogic,1,'first');
leftIdx = find(r <= eps & ~surfLogic);
surfEnd = find(surfLogic,1,'last');
switch condition
    case 'spread stick'
        dof = 2*length(z);
        gVec = zeros(dof,1);
        gLogic = zeros(dof,1);
        gLogic(2*leftIdx - 1) = 1;
        gLogic(2*startSurf-1:2*surfEnd) = 1;
        hVec = ones(dof,1);
        hVec(2*startSurf-1:end) = 0;
        hVec(2*leftIdx) = 1;
    case 'spread slip'
        dof = 2*length(z);
        gVec = zeros(dof,1);
        gLogic = zeros(dof,1);
        gLogic(2*startSurf:2:end) = 1;
        gLogic(1) = 1;
        hVec = ones(dof,1);
        hVec(2*startSurf:2:end) = 0;
    case 'top plate'
        dof = 2*length(z);
%         topIdx =  find(z>=(max(z)-eps),1,'last');
        gVec = zeros(dof,1);
        gLogic = zeros(dof,1);
        gLogic(2*startSurf:2:end) = 1;
        gLogic(1) = 1;
        hVec = ones(dof,1);
        hVec(2*startSurf:2:end) = 0;
    case 'pressure'
        dof = length(z);
        gVec = zeros(dof,1);
        gLogic = zeros(dof,1);
        hVec = ones(dof,1);
end
gCell = {gVec,gLogic};
hCell = {hVec};

end