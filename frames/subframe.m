function x = subframe(x,y)
    for m = 1:numel(x)
        for n = 1:numel(x{m})
            x{m}{n} = x{m}{n} - y{m}{n};
        end
    end
end