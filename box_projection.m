function x=box_projection(x,lower_bound,upper_bound)
% Projection of 2D data x onto its constraints defined by lower and upper bounds
% For projections on quantization constraint sets
x = min(cat(3, max(cat(3,lower_bound,x),[],3), upper_bound),[],3);
end