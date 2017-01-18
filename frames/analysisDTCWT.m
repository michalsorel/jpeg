function pyr = analysisDTCWT(x,Nsc,par2,par3)
%
[Yl,Yh] = dtwavexfm2(x,Nsc,par2,par3);
pyr = cell(2*Nsc+1,1); 
for i = 1:numel(Yh) % need to decompose to complex and real parts
    pyr{2*i-1} = real(Yh{i});
    pyr{2*i} = imag(Yh{i});
end
pyr{2*Nsc+1} = Yl;