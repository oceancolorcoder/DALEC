function mw=Index_w(wv, T, S)
% Calculate the refractive index of water
% wv: wavelength in nanometer
% T: temperature in degree
% S: Salinity in per thousand
% mw(wv,T,S)=n0+(n1+n2T+n3T^2)S+n4T2+(n5+n6S+n7T)/wv+n8/wv^2+n9/wv^3;

n0=1.31405;
n1=1.779e-4;
n2=-1.05e-6;
n3=1.6e-8;
n4=-2.02e-6;
n5=15.868;
n6=0.01155;
n7=-0.00423;
n8=-4382;
n9=1.1455e6;

n0_4=n0+(n1+n2*T+n3*T^2).*S+n4*T^2;
n5_7=n5+n6*S+n7*T;

mw=n0_4+n5_7*(wv.^-1)+n8*(wv.^-2)+n9*(wv.^-3);
end

