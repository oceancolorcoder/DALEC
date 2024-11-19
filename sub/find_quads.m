function loc = find_quads(quads,zen,azm)
% find locations in quads
zen = zen(:);
azm = azm(:);
if ~isequal(numel(zen),numel(azm))
    if isscalar(zen)
        zen = zen*ones(size(azm));
    elseif isscalar(azm)
        azm = azm*ones(size(zen));
    else
        error('Non-matching sizes');
    end
end

loc = zen;

for i=1:numel(zen)
    tmp = ((quads.zen-zen(i)).^2+(quads.azm-azm(i)).^2).^0.5;
    [~,loc(i)]=min(tmp);
end
end