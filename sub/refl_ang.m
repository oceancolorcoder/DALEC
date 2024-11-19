function [theta_r,phi_r,ang]=refl_ang(theta0,phi0,theta_i,phi_i)
% Calculate the fresnel reflection angle given the facet normal
% [theta_r,phi_r,ang]=refl_ang(theta0,phi0,theta_i,phi_i,m)
% theta_r and phi_r are zenith and azimuth angles of the reflected ray
% theta0 and phi0 are zenith and azimuth angles of the normal of facet
% theta_i and phi_i are zenith and azimuth angles of the incident ray
% ang: the incident angle

% Reference: Light and Water: Radiative Transfer in Natural Waters
%            Mobley, 1994 (Page 156, Fig. 4.2)

% Vector expression of incident and normal
% inc=[sin(theta_i).*cos(phi_i) sin(theta_i).*sin(phi_i) cos(theta_i)]';
% nor=[sin(theta0).*cos(phi0) sin(theta0).*sin(phi0) cos(theta0)]';
inc = my_sph2cart(phi_i, theta_i, 1);
nor = my_sph2cart(phi0, theta0, 1);
ref=zeros(size(inc));

% angle between incident and normal
cosw = sum(bsxfun(@times,nor,inc),2);
% cosw=sum(inc.*nor,2);
ang = acos(cosw);
ind = ang>pi/2;
ang(ind) = pi - ang(ind);
% weight=Fresnel(m,acos(cosw)');

% ref = inc - 2*cosw*nor
ref = bsxfun(@minus,inc,2*bsxfun(@times,nor, cosw));

[phi_r,theta_r,~] = my_cart2sph(ref);
% theta_r = acos(ref(3,:))';
% phi_r = atan(ref(2,:)./ref(1,:))';
% ind = find(ref(1,:)==0);
% if ref(2,ind)>=0
%    phi_r(ind) = pi/2;
% else
%    phi_r(ind) = 3/2*pi;
% end
% ind = find(ref(1,:)<0);
% phi_r(ind) = pi + phi_r(ind);
end