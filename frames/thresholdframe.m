function [res, u_v] = thresholdframe(pyr,mutau,threshtype,threshf)
% threshold wavelet coefficients stored in a cell array of cell arrays of matrices
    res = cell(size(pyr));
    u_v = zeros(size(pyr));
    for ii = 1:numel(pyr)   
        for iii = 1:numel(pyr{ii})
            res{ii}{iii} = thresholding(pyr{ii}{iii},mutau,threshtype,threshf);                        
        end
    end
    if nargout > 1
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
