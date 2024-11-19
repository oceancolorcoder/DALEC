function p = normalcdf(x0, u,sigma)
% For normal distribution, the CDF is
% cdf(x) = 1/2*(1+erf((x-u)/sigma/sqrt(2)))
% p = integral(@(x)normaldp(x,u,sigma),-inf,x0);
p = 0.5*(1+erf((x0-u)/sigma/2^0.5));
end