function vec = gen_vec(zens,azms)
% generate vectors from permutation of zenith and azimuth angles
[zens,azms] = ndgrid(zens,azms);
zens = zens(:);
azms = azms(:);
% vector expression
vec = my_sph2cart(azms,zens,1);
end