function x = synthesisDWT(pyr,par2)
%
% Image reconstruction form pyramidal decomposition using
% Matlab Wavelet Toolbox
% par2 = 'name_of_wavelet' e.g.'db1' see help of wfilters.m

% Reorganize coefficients into Matlab Wavelet Toolbox format
size_x=size(pyr{1});
pyr=pyr(2:end);
% S=zeros(length(pyr)+1,2);
S=zeros(length(pyr),2);

for n=1:length(pyr)
    S(n+1,:)=[size(pyr{n},1) size(pyr{n},2)];
    pyr{n}=reshape(pyr{n},1,[]);
end
C=[pyr{end:-1:1}];
S=flipud(S);
% S(n+1,:)=S(end-1,:)*2;
S(n+1,:)=size_x;

% inverse DWT (Matlab Wavelet Toolbox)
x = waverec2(C,S,par2);