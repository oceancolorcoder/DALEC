function vec = gen_vec_quad(zen,du,azm,dphi,num)
% generate vectors for non-polar quads
narginchk(4,5)
if nargin <5
    num = 10;
end
half_azm = linspace(-dphi/2,dphi/2,num);
half_zen = linspace(-du/2/sin(zen), du/2/sin(zen),num);
vec = gen_vec(zen+half_zen, azm+half_azm);
end