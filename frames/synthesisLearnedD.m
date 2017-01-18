function x = synthesisLearnedD(pyr,dict)

row = size(pyr{1},1);
col = size(pyr{1},2);
% scalar = diag(dict'*dict);
scalar = squeeze(sum(sum(dict.*dict,1),2));
filterSize = size(dict,1);
no_of_filters=size(dict, 3);

x=filter2(rot90(dict(:,:,1),2), pyr{end}/scalar(1), 'full');

for n = 2:no_of_filters
    kernel = rot90(dict(:,:,n),2);
    x = x + filter2(kernel, pyr{1}(:,:,n-1)/scalar(n-1), 'full');
end

%% reconstruction
row=row+filterSize-1;
col=col+filterSize-1;
for k = 1:filterSize
    for k2 = 1:filterSize
        if (k == 1)&(k2 == 1)
            mmask = ones(row, col);
        else
            if k == 1
                temp = zeros(row, col);
                temp(:, k2:col-filterSize+k2-1) = 1;
                mmask = mmask + temp;
            elseif k2 == 1
                temp = zeros(row,col);
                temp(k:row-filterSize+k-1, :)=1;
                mmask = mmask + temp;
            else
                temp = zeros(row, col);
                temp(k:row-filterSize+k-1, k2:col-filterSize+k2-1) = 1;
                mmask = mmask + temp;                
            end
        end
    end
end
mmask  = double(mmask);
x = x./mmask;