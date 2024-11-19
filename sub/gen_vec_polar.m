function vec = gen_vec_polar(zen,sun05,num)
% generate vectors for the polar cap quad and sun disk
% By convention, the sun disk is at XZ plane, i.e., azimuth = 0.
% zen: the zenith angle of Sun

narginchk(2,3)
if nargin <3
    num = 10;
end

%
phi = linspace(0,2*pi,num);
tmp = [sin(sun05)*cos(phi); sin(sun05)*sin(phi); ones(size(phi))*cos(sun05)];
tmp = [[0 0 1]' tmp]; % add the center vector
% tmp = tmp/sqrt(1+tan(sun05)^2);
% rotation about Y axis
Ry = [cos(zen),0,sin(zen);0,1,0;-sin(zen),0,cos(zen)];
vec = (Ry*tmp)';
%}

end