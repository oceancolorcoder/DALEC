function [azm, zen, r] = my_cart2sph(n)
[azm,zen,r] = cart2sph(n(:,1),n(:,2),n(:,3));
zen = pi/2 - zen;
end