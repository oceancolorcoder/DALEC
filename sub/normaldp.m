function p = normaldp(x, u,sigma)
p = (1/(sqrt(2*pi)*sigma))*exp(-((x-u).^2)./(2*sigma^2));
end