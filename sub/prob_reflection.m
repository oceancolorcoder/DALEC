function [prob,ang] = prob_reflection(inc,refl,wind)
% estimate the probablity of facets that would reflect incident light into
% the direction defined by refl under winds
% inc:      incident light vector
% refl:     reflected light vector
% wind:     wind speeds
% prob:     probability
% ang:     reflection angle

% the norm vector of the facet that would reflect sun light into sensor
n = bsxfun(@minus, refl, inc);
% n = bsxfun(@rdivide, n, sum(abs(n).^2,2).^(1/2));
n = bsxfun(@rdivide, n, vec_length(n));

% the zenith and azimuth angles of the facets
[azm_n,zen_n,r] = my_cart2sph(n);

% convert facet zenith angle to slopes
slope = tan(zen_n);

% estimate wind-roughned probability of facets
% sigma2 = 0.003 + 0.00512*wind;
% sigma = sigma2^0.5;
sigma = cox_munk(wind);
% p1 = normcdf(max(slope),0,sigma) - normcdf(min(slope),0,sigma);
% !!! see document On the Cox and Munk
%
sigma = sigma/sqrt(2);
p1 = rayleighcdf(max(slope),sigma)-rayleighcdf(min(slope),sigma);
%} !!! 
% azimuth angle ranges from -180 to 180. Need to treat the cases when the
% azimuth angles cover both positive ang negative ranges.
% case 1: -2 -1 1 2
% case 2: -179, -178, 178, 179
% case 3: -179 -120 2 5 130 178
% cases 1 and 2: the range should be 4
% case 3: the range should be 357
azm_nx = max(azm_n);
azm_nn = min(azm_n);
if azm_nx*azm_nn >0 % not an issue
    p2 = (azm_nx-azm_nn)/2/pi;
elseif any(abs(azm_n)<pi/2) % cases 1 and 3 
    p2 = (azm_nx-azm_nn)/2/pi;
else % case 2
    ind = azm_n<0;
    azm_n(ind) = azm_n(ind)+2*pi;
    azm_nx = max(azm_n);
    azm_nn = min(azm_n);
    p2 = (azm_nx-azm_nn)/2/pi;
end
prob = 2*p1*p2; % factor 2 accounts for 180 degree ambiguity

% incident angle
cosw = sum(bsxfun(@times,n,refl),2);
ang = acos(cosw);
ind = ang>pi/2;
ang(ind) = pi - ang(ind);
ang = mean(ang);
end