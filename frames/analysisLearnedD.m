function pyr = analysisLearnedD(x,dict)
%

[row,col] = size(x);
filterSize = [size(dict,1) size(dict,2)];
no_of_filters=size(dict, 3);
pyr={zeros(row-filterSize(1)+1,col-filterSize(2)+1,no_of_filters-1),...
    zeros(row-filterSize(1)+1,col-filterSize(2)+1)};

pyr{2} = filter2(dict(:,:,1), x, 'valid'); % First is low pass?
for n = 2:no_of_filters;
    kernel=dict(:,:,n);
    pyr{1}(:,:,n-1) = filter2(kernel, x, 'valid');
end