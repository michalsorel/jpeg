function x = synthesisDTCWT(pyr,par2,par3)
%
Yl = pyr{end}; % last one is always low-pass
hs = (numel(pyr)-1)/2;
Yh = cell(hs,1);
for ii = 1:hs
    Yh{ii} = pyr{2*ii-1} + i*pyr{2*ii};
end
x = dtwaveifm2(Yl,Yh,par2,par3); 



