function ref = sw_Fresnel(wv,ang,T,S)
% Calcualte Fresnel reflectance for seawater
% wv: wavelength in nm
% ang: reflectanc angle
% T: temperature in Celsius
% S: salinity in PSU
m = Index_w(wv,T,S);
ref = Fresnel(m,ang);
end