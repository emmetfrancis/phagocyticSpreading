function adhForce = adhForceDiscreteAxisymm_r4(r,z,maxDist,r0Adh,f0,NLoc,rEdge, TForce)

adhForce = zeros(2,length(z));
spacing = mean(diff(NLoc));
% rEdge = r(end);
% figure
for i = 1:length(z)
    if abs(z(i)-maxDist) < 1e-12 % no region to integrate over
        adhForce(:,i) = [0, 0];
    elseif abs(z(i)) < 1e-12 % consider in contact
        adhForce(:,i) = [0,0];
    else % need to integrate over region
        r0 = r0Adh;% * TForce(i)^.1;
        rm = r(i);
        zm = z(i);
        for j = find(NLoc > rEdge)
            NLocCur = NLoc(j);
            DCur = sqrt((rm-NLocCur)^2 + zm^2);
            fCur = f0*r0^4*(1/DCur^4 - r0^3 * (1/DCur^7));
            adhForce(1,i) = adhForce(1,i) + (1/spacing)*(NLocCur-rm)*NLocCur*fCur/DCur;
            adhForce(2,i) = adhForce(2,i) + (1/spacing)*(-zm)*NLocCur*fCur/DCur;
        end
    end
end