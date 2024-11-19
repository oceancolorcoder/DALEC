function n = my_sph2cart(azm,zen,r)
[x,y,z] = sph2cart(azm,pi/2-zen,r);
n = [x(:) y(:) z(:)];
end