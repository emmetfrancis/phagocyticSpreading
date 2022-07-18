function [mCurv,pCurv,phi] = polyCurv(s,r,z,polyOrder,numFitIndices)
% polyOrder must be at least 2 (because second derivative)

if polyOrder < 2
    polyOrder = 2;
end
if isempty(s)
    s = zeros(1,length(r));
    for i = 2:length(r)
        s(i) = s(i-1) + sqrt((r(i)-r(i-1))^2 + (z(i)-z(i-1))^2);
    end
end
% create matrices to save coefficients for polynomial fits
zp = zeros(length(z), polyOrder+1);
dzp = zeros(length(z), polyOrder);
d2zp = zeros(length(z), polyOrder-1);
rp = zeros(length(r), polyOrder+1);
drp = zeros(length(r), polyOrder);
d2rp = zeros(length(r), polyOrder-1);
halfInt = floor(numFitIndices/2);
for i = 1:length(r)
    if i < halfInt+1 % handle points at start
        if i >= polyOrder+1
            fitIndices = 1:i+halfInt;
        else
            continue
        end
    elseif i > length(r)-halfInt % handle points at end
        if (length(r)-i+1) >= polyOrder+1
            fitIndices = i-halfInt:length(r);
        else
            continue
        end
    else
        fitIndices = (i-halfInt):(i+halfInt);
    end
    if length(fitIndices) < polyOrder + 1
        continue
    end
    zp(i,:) = polyfit(s(fitIndices), z(fitIndices), polyOrder);
    dzp(i,:) = polyder(zp(i,:));
    d2zp(i,:) = polyder(dzp(i,:));
    rp(i,:) = polyfit(s(fitIndices), r(fitIndices), polyOrder);
    drp(i,:) = polyder(rp(i,:));
    d2rp(i,:) = polyder(drp(i,:));
end
normal = zeros(2, length(r));
phi = zeros(1, length(r));
mCurv = zeros(1, length(phi));
for i = 1:length(r)
    dzCur = polyval(dzp(i,:), s(i)); %dz/ds
    d2zCur = polyval(d2zp(i,:), s(i));
    drCur = polyval(drp(i,:), s(i)); %dr/ds
    d2rCur = polyval(d2rp(i,:), s(i));
    normal(:,i) = [drCur; -dzCur] ./ (sqrt(drCur^2+dzCur^2)); %unit normal
    mCurv(i) = -drCur*d2zCur + dzCur*d2rCur;
end
phi = atan2(normal(2,:), normal(1,:));
% phi(phi < -pi/2) = phi(phi < -pi/2) + 2*pi;

% values at end and start set to value of nearest neighbor that can be
% computed (could try some sort of extrapolation)
firstComp = find(~isnan(phi),1,'first');
lastComp = find(~isnan(phi),1,'last');
for i = 1:halfInt
    if isnan(phi(i))
        phi(i) = phi(firstComp);
        mCurv(i) = mCurv(firstComp);
    end
end
for i = length(r)-halfInt+1:length(r)
    if isnan(phi(i))
        phi(i) = phi(lastComp);
        mCurv(i) = mCurv(lastComp);
    end
end
mCurv(1:halfInt) = mCurv(halfInt+1);
mCurv(length(r)-halfInt+1:end) = mCurv(length(r)-halfInt);
phi(1:halfInt) = phi(halfInt+1);
phi(length(r)-halfInt+1:end) = phi(length(r)-halfInt);

% smooth mCurv
% mCurv = smooth(s,mCurv,9);

pCurv = sin(phi) ./ r; % by definition
pCurv(1:halfInt) = pCurv(halfInt+1);
pCurv(length(r)-halfInt+1:end) = pCurv(length(r)-halfInt);
% pCurv(r<.001*max(r)) = mCurv(r<.001*max(r));

end