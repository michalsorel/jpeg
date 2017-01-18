function r = snr(original,degraded)
%SNR using original and degraded images
%
%   function r = snr(original,degraded)
%
%Remark: degraded can be cell array, then r is cell array
%
%Example: snr(q,noise(q,20)) results in cca 20.
%
%Michal Sorel (c) 2005, 2009

original=double(original);
degraded=double(degraded);
if ~iscell(degraded)
    r = 10*log10(var(original(:))...
                  / var(original(:)-degraded(:)));
else
    r = cell(size(degraded));
    for i = 1:length(degraded)        
        r{i} = 10*log10(var(original{i}(:))...
                  / var(original(:)-degraded{i}(:)));
    end
end