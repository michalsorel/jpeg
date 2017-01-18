function [res, u_v] = threshhighpass(pyr,mutau,threshtype,threshf)
% threshold high-passes of the pyramid
    res = cell(size(pyr));
    u_v = zeros(size(pyr));
    for ii = 1:numel(pyr)-1 % without the last low-pass band                
        res{ii} = thresholding(pyr{ii},mutau,threshtype,threshf);
        u_v(ii)=(norm(res{ii}(:)-pyr{ii}(:)))^2;
    end
    res{numel(pyr)} = pyr{numel(pyr)};    
    % Low-pass only
    res{ii+1} = thresholding(pyr{ii+1},mutau,threshtype,threshf);
    u_v(ii+1)=(norm(res{ii+1}(:)-pyr{ii+1}(:)))^2;
    u_v = sqrt(u_v);
end
