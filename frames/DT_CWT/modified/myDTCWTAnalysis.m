% Dual-Tree Complex Wavelet Transform (Analysis)
% Adapted from Kingsbury's original code to match Simoncelli's Matlab
% Pyrtools
% Author: Luis Mancera
% Date: 21/02/2006
% Revised: Javier Portilla in March/April 2008
%
% [pyr,pind] = myDTCWTAnalysis(X, par1, par2, par3)
% 
%   Inputs:
%       X  - Image
%       par1 - Number of levels
%       par2 - biort:
%               - 'antonini' => Antonini 9.7 tap filters
%               - 'legall'     => LeGall 5,3 tap filters.
%               - 'near_sym_a' => Near-Symmetric 5,7 tap filters.
%               - 'near_sym_b' => Near-Symmetric 13,19 tap filters.
%       par3: qshift: 
%               - 'qshift_06' => Quarter Sample Shift Orthogonal 
%                           (Q-Shift) 10,10 tap filters, 
%                          (only 6,6 non-zero taps).
%               - 'qshift_a' =>  Q-shift 10,10 tap filters,
%                          (with 10,10 non-zero taps, unlike qshift_06).
%               - 'qshift_b' => Q-Shift 14,14 tap filters.
%               - 'qshift_c' => Q-Shift 16,16 tap filters.
%               - 'qshift_d' => Q-Shift 18,18 tap filters.
%
function [pyr,pind] = myDTCWTAnalysis(X, par1, par2, par3)
   
    Nor = 6;
    [Yl,Yh] = dtwavexfm2(X,par1,par2,par3);
        
    % Insert the subbands into pyr
    pyr = [];
    pind = [];
    for s=1:length(Yh)
        for w=1:Nor,
            aux = Yh{s}(:,:,w);
            pyr = [pyr; aux(:)];
            pind = [pind; size(aux)];
        end
    end
    
    % Insert the low pass residual
    pyr = [pyr; Yl(:)];
    pind = [pind; size(Yl)];
   
    % Finally, separate real from imaginary part
    pyr = [real(pyr); imag(pyr)];
    pind = [pind; pind];
return
