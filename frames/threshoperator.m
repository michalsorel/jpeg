function [x_r,varargout] = threshoperator(z,analysis,synthesis,mutau,par)
% [x_r,{A, B}] = threshoperator(z,analysis,synthesis,mutau,par)
% A=||pyr_thresh-pyr||^2, B=(||pyr||_p)^p

size_z=size(z);
x_r = zeros(size_z,class(z));
no_of_frames = numel(analysis);
value=0;
norm_a2=0;
pyr_out=cell(1,no_of_frames);
for m = 1:no_of_frames
    pyr = analysis{m}(z);
    pyr = mulcell(pyr,sqrt(no_of_frames)/no_of_frames);
    new = threshhighpass(pyr,mutau,par.threshtype{m},par.threshf{m});
    synth=synthesis{m}(new);
%     x_r = x_r + synthesis{m}(new);
    x_r = x_r + synth(1:size_z(1),1:size_z(2));

    if (nargout-1)~=0
        for n=1:length(pyr)
            resid=pyr{n}-new{n};
            value=value+sum(resid(:).^2);
            norm_a2=norm_a2 + sum(abs(new{n}(:)).^par.threshtype{m});
        end
        pyr_out{m}=new;
        var_out=[value norm_a2 {pyr_out}];
    end
end
x_r = sqrt(no_of_frames)/no_of_frames * x_r;

for k=1:nargout-1
    varargout(k) = var_out(k);
end

end