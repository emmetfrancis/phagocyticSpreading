function [SA,vol] = SAVolCalc(r,z,method);
% calculate SA and volume for an axisymmetric body

sParam(1) = 0;
for i = 2:length(r)
    sParam(i) = sParam(i-1) + sqrt((r(i)-r(i-1)).^2 + (z(i)-z(i-1)).^2);
end
switch method
    case 'trapz'
        SA = trapz(sParam,2*pi*r);
        vol = -trapz(z,pi*r.^2); %because z goes from high to low, flip integral
    otherwise
        error('Other SA and volume integration methods not coded yet')
        SA = NaN;
        vol = NaN;
end

end