function y = rayleighcdf(x,s)
% Cumulative distribution function for Rayleigh distribution
t = (x./s).^2;
y = 1-exp(-t/2);
end