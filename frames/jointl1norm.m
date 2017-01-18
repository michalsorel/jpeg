function r = jointl1norm(x)
    r = 0;
    for m = 1:numel(x{1})
        for n = 1:numel(x{1}{m})
            nrm = 0;
            for gi = 1:numel(x)
                nrm = nrm + x{gi}{m}{n}.^2;
            end
            r = r + sum(sqrt(nrm(:)));            
        end
    end
end
