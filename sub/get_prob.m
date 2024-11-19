function [prob,angr_sky] = get_prob(wind,vec,quads)
% prob: probability of sky light reflected into the sensor
% angr_sky: reflection angle
prob = nan(numel(quads.zen),numel(wind)*size(vec,1));
angr_sky = prob;
k = 0;
for i=1:numel(wind)
    for j=1:size(vec,1)
        k = k+1;
        [prob(:,k),angr_sky(:,k)] = sky_light_reflection2(wind(i),vec(j,:),quads);
    end
end
end