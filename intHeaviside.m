function [int_h] = intHeaviside(t)
%INTHEAVISIDE Integral of the Heaviside step: returns t for t>=1, 0 otherwise
int_h = zeros(size(t));
idx = t >= 1;
int_h(idx) = t(idx);
end