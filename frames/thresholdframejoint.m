function [res, u_v] = thresholdframejoint(pyr,mu2tau)
% threshold wavelet coefficients stored in a cell array of cell arrays of
% matrices and grouped into an additional cell level 
% All frames must have the same dimensions on all levels
    res = cell(size(pyr));
    u_v = zeros(size(pyr));
    grouplen = numel(pyr);
    for ii = 1:numel(pyr{1})   
        for iii = 1:numel(pyr{1}{ii})
            nrm = 0;
            for gi = 1:grouplen
                nrm = nrm + pyr{gi}{ii}{iii}.^2;
            end
            for gi = 1:grouplen
                res{gi}{ii}{iii} = soft_thresh_group(pyr{gi}{ii}{iii},sqrt(nrm),1/mu2tau(1));                        
            end
        end
    end
    if nargout > 1
        error('Not implemented');
        for ii = 1:numel(pyr)   
               u_v(ii) = 0;
            for iii = 1:numel(pyr{ii})
                u_v(ii)=u_v(ii) + (norm(res{ii}{iii}(:)-pyr{ii}{iii}(:)))^2;
            end
        end
            u_v = sqrt(u_v);
    end
    %res{numel(pyr)} = pyr{numel(pyr)};        
end

function f_out = soft_thresh_group(f,normf,theta) % group soft thresholding function    
    %f_out = sign(f).*max(abs(f)-theta,0);
    ind = normf > 0;
    f_out = zeros(size(f));
    f_out(ind) = f(ind)./normf(ind).*max(normf(ind)-theta,0);
    f_out(~ind)=0;
end