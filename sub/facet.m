function [theta,phi]=facet(wind,len)
% Calculate the zenith and azimuth angle of sea surface normal given wind speed
% [theta,phi]=Pfacet(wind)
% wind: wind speed
% The slope (= tan(theta)) obeys normal distribution, whose variance sigma^2
% depends on wind speed through
% sigma^2 = 0.003 + 0.00512*wind

sigma = cox_munk(wind);
Sslope = randn(len,1)*sigma; % sample of real slope
theta = pi/2+abs(atan(Sslope));
phi = 2*pi*rand(len,1);
end