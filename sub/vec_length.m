function al = vec_length(a)
% the length of vector a
al = sum(abs(a).^2,2).^(1/2);
end