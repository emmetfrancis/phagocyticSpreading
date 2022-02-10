function [mCurv,pCurv,phi] = localCurv(r,z,polyOrder,numFitIndices)

% polyOrder must be at least 2 (because second derivative)

if polyOrder < 2
    polyOrder = 2;
end
if numFitIndices < polyOrder + 1
    numFitIndices = polyOrder + 1;
end
s = zeros(1,length(r));
for i = 2:length(r)
    s(i) = s(i-1) + sqrt((r(i)-r(i-1))^2 + (z(i)-z(i-1))^2);
end
% create matrices to save coefficients for polynomial fits
zp = zeros(length(z), polyOrder+1);
dzp = zeros(length(z), polyOrder);
d2zp = zeros(length(z), polyOrder-1);
rp = zeros(length(r), polyOrder+1);
drp = zeros(length(r), polyOrder);
d2rp = zeros(length(r), polyOrder-1);
halfInt = floor(numFitIndices/2);
dsVals = zeros(1,length(r));
dsVals(1:end-1) = smooth(diff(s),numFitIndices);
for i = 1:length(r)-1
    sStart = s(i) - halfInt*dsVals(i);
    sEnd = s(i) + halfInt*dsVals(i);
    if sStart < 0 % handle points at start
        if sEnd >= polyOrder*dsVals(i)
            sStart = 0;
            sEnd = floor(sEnd/dsVals(i))*dsVals(i);
        else
            continue
        end
    elseif sEnd > s(end) % handle points at end
        if sStart <= s(end)-polyOrder*dsVals(i)
            sStart = s(end) - floor((s(end)-sStart)/dsVals(i))*dsVals(i);
            sEnd = s(end);
        else
            continue
        end
    end
    sFit = sStart:dsVals(i):sEnd;
    zFit = interp1(s,z,sFit);
    rFit = interp1(s,r,sFit);
    zp(i,:) = polyfit(sFit, zFit, polyOrder);
    zDer = polyder(zp(i,:));
    dzp(i,end-length(zDer)+1:end) = zDer; % in case first coeff is zero (otherwise throws error)
    z2Der = polyder(dzp(i,:));
    d2zp(i,end-length(z2Der)+1:end) = z2Der;
    rp(i,:) = polyfit(sFit, rFit, polyOrder);
    rDer = polyder(rp(i,:));
    drp(i,end-length(rDer)+1:end) = rDer;
    r2Der = polyder(drp(i,:));
    d2rp(i,end-length(r2Der)+1:end) = r2Der;
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
    if isnan(phi(i)) || isinf(phi(i))
        phi(i) = phi(firstComp);
        mCurv(i) = mCurv(firstComp);
    end
end
for i = length(r)-halfInt+1:length(r)
    if isnan(phi(i)) || isinf(phi(i))
        phi(i) = phi(lastComp);
        mCurv(i) = mCurv(lastComp);
    end
end
% mCurv(1:halfInt) = mCurv(halfInt+1);
% mCurv(length(r)-halfInt+1:end) = mCurv(length(r)-halfInt);
% phi(1:halfInt) = phi(halfInt+1);
% phi(length(r)-halfInt+1:end) = phi(length(r)-halfInt);

% smooth mCurv
% mCurv = smooth(s,mCurv,9);

pCurv = sin(phi) ./ r; % by definition
% pCurv(1:halfInt) = pCurv(halfInt+1);
% pCurv(length(r)-halfInt+1:end) = pCurv(length(r)-halfInt);
% pCurv(r<.001*max(r)) = mCurv(r<.001*max(r));
if any(isnan(phi))
    phi(isnan(phi)) = pi;
end

end