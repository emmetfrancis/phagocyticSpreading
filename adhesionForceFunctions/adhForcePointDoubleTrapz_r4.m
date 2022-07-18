function adhForce = adhForcePointDoubleTrapz_r4(r,z,maxDist,r0Adh,f0,N,rEdge,repulseVal)

adhForce = zeros(2,length(z));
% rEdge = r(end);
% figure
for i = 1:length(z)
    if abs(z(i)) < 1e-12 % consider in contact
        adhForce(:,i) = [0,0];
    else%if z(i) > r0Adh*((7/4)^(1/3)) % need to integrate over region
        r0 = r0Adh;% * TForce(i)^.1;
        rm = r(i);
        zm = z(i);
%         rSurf = rEdge:.001:rEdge+100*r0;
        rSurf = [rEdge,rEdge+logspace(-6,log10(50*r0),50)];
        refAngle = 50*r0/rm;
        if abs(refAngle) > pi
            refAngle = pi;
        end
        angSurf = logspace(log10(pi/10000),log10(refAngle/2),50);
        angSurf = [-angSurf(end:-1:1),0,angSurf];
%         angSurf = -refAngle/2:pi/1000:refAngle/2;
        [rMesh,angMesh] = meshgrid(rSurf,angSurf);
        DMesh = sqrt((rMesh.*cos(angMesh)-rm).^2 + (rMesh.*sin(angMesh)).^2 + zm.^2);
        FrMesh = (rMesh.*cos(angMesh)-rm) .* (N*(1./DMesh).^5 - r0^3 * repulseVal * (1./DMesh).^8) .* rMesh;
        FzMesh = -zm .* (N*(1./DMesh).^5 - r0^3 * repulseVal * (1./DMesh).^8) .* rMesh;
%         FrMesh = (rMesh.*cos(angMesh)-rm) .* (N*(1./DMesh).^4 - r0^3 * repulseVal * (1./DMesh).^7) .* rMesh;
%         FzMesh = -zm .* (N*(1./DMesh).^4 - r0^3 * repulseVal * (1./DMesh).^7) .* rMesh;
        adhForce(1,i) = (f0*r0^4)*trapz(angSurf,trapz(rSurf,FrMesh,2));
        adhForce(2,i) = (f0*r0^4)*trapz(angSurf,trapz(rSurf,FzMesh,2));
    end
end

end