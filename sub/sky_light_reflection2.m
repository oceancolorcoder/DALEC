function [prob,ang] = sky_light_reflection2(wind,sensor,quads)

% wind: wind speed (m/s);
% sensor: The vector of reflected light to be measured by the sensor
% 
% quads: the quads of sky light
%
% Xiaodong Zhang
% 2015-05-18:
% 2015-05-30:   use probability instead of Monte Carlo


% initialize
prob = quads.zen;
ang = prob;

% polar quad, 1st in quads
zen0 = quads.zen(1);
% generate sky vector
num = 100;
p_vec = gen_vec_polar(zen0,quads.sun05,num);
% -p_vec represent vectors coming from the sky
[prob(1),ang(1)] = prob_reflection(-p_vec,sensor,wind);

% non-polar quads
num = 10; % the number of individual vectors
%{
half_azm = linspace(-quads.dphi/2,quads.dphi/2,num);
for i=2:length(prob)
    half_zen = linspace(-quads.du/2./sin(quads.zen(i)), quads.du/2./sin(quads.zen(i)),num);
    sky = gen_vec(quads.zen(i)+half_zen, quads.azm(i)+half_azm);
    [prob(i),ang(i)] = prob_reflection(sky,sensor,wind);
end
%}
for i=2:length(prob)
    sky = gen_vec_quad(quads.zen(i),quads.du,quads.azm(i),quads.dphi,num);
    [prob(i),ang(i)] = prob_reflection(-sky,sensor,wind);
end
end

