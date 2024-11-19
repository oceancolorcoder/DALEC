function [R,R12,R33]=Fresnel(m,ang)
% This function calculate the Fresnel reflectances for electric vector
% parallel (Rp), perpendicular (Rr) and unpolarized incident light.
% ang is incident angle.
% m is relative index of refrax.
% 
% The reflection matrix = 
% [R11, R12, 0; R12, R11, 0; 0, 0, R33]
% Only accounts for I, Q, U and ignore the V component.
%% Revision History
% 2016-07-10:   1st version, just compute R11, i.e, R
% 2016-12-14:   add other reflection matrix elements R12 and R33
%               Also found an error in the previous equaiton for Rp1

%%
ang = ang(:); % column vector
m = m(:)'; % row vector

cosang=abs(cos(ang)); % cosine of incident angle
sinangr=sin(ang)*(1./m); % sine of refraction angle
cosangr=(1-sinangr.^2).^0.5; % cosine of refraction angle

% reflection coefficient for perpendicular incident light
tmp = bsxfun(@times,cosangr,m);
Rr1 = bsxfun(@minus,cosang,tmp)./bsxfun(@plus,cosang,tmp);
% Rr1=(cosang-m*cosangr)./(cosang+m*cosangr);

% reflection coefficient for parallel incident light
tmp = bsxfun(@times,cosang,m);
% this was previous one
% Rp1 = bsxfun(@minus,cosangr,tmp)./bsxfun(@plus,cosangr,tmp);    
Rp1 = bsxfun(@minus,tmp, cosangr)./bsxfun(@plus,cosangr,tmp);    
% Rp1=(cosangr-m*cosang)./(cosangr+m*cosang);

Rr=abs(Rr1).^2; % reflectance for perpendicular incident light

Rp=abs(Rp1).^2; % reflectance for parallel incident light

R=(Rr+Rp)/2;
R12=(Rp-Rr)/2;
R33=real(Rr1.*conj(Rp1));
end

