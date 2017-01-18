function sum_x=sum_L1_frame(x,analysis)
if ~iscell(x)
    if size(x,3)>1
        x=squeeze(mat2cell(x,size(x,1),size(x,2),ones(1,size(x,3))));
    else
        x=squeeze(mat2cell(x,size(x,1),size(x,2)));
    end
end

f=@(a)sum(abs(a(:)));
L1pyr=zeros(length(x),length(analysis));
for p=1:length(x) % over image components
    for m = 1:length(analysis) % over different frames (DT-CWT,Haar)
        pyr = analysis{m}(x{p});
        L1pyr(p,m)=sum(cellfun(f,pyr));
    end
end
sum_x=sum(L1pyr(:));