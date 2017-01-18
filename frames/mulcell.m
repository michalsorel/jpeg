function a = mulcell(a,b)
% multiplication of cell arrays

if iscell(b)
    for i = 1:numel(a)
        a{i} = a{i}*b{i};
    end
else
    for i = 1:numel(a)
        a{i} = a{i}*b;
    end    
end