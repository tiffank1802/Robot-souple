function [h] = Heaviside(t)
%HEAVISIDE Vectorized Heaviside step function: 0 for t<1, 1 for t>=1
h = zeros(size(t));
% Use linear indexing so this works for row/column/matrix inputs
idx = t >= 1;
h(idx) = 1;
end