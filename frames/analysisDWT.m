function pyr = analysisDWT(x,Nsc,par2)
%
% Image wavelet decomposition using Matlab Wavelet Toolbox
% Nsc = number of levels
% par2 = 'name_of_wavelet' e.g.'db1' see help of wfilters.m

% DWT using Matlab Wavelet Toolbox
[C,S] = wavedec2(x,Nsc,par2);

% Reorganize coefficients into our format
pyr = cell(Nsc+1,1);
st_ind=S(1,1)*S(2,2);
pyr{Nsc+1}=reshape(C(1:st_ind),S(1,1),S(1,2));

for n=2:size(S,1)-1
    n_el=S(n,1)*S(n,2)*3;
    pyr{-n+size(S,1)}=reshape(C((st_ind+1):(st_ind+n_el)),S(n,1),S(n,2),3);
    st_ind=st_ind+n_el;
end
pyr=[{zeros(size(x,1),size(x,2))}; pyr];
end