% Dual-Tree Complex Wavelet Transform (Synthesis)
% Adapted from Kingsbury's original code to match Simoncelli's Matlab
% Pyrtools
% Author: Luis Mancera
% Date: 21/02/2006
% Revised: Javier Portilla in March/April 2008
%
% X = myDTCWTSynthesis(pyr, pind, par2, par3)
%
%   Inputs:
%       pyr  - Dual-Tree Complex Wavelet Transform. Composed by
%             Yl -> The real lowpass image from the final level
%             Yh -> A cell array containing the 6 complex highpass 
%                   subimages for each level.
%       pind    - pyramid indices
%       par2 - biort:
%               - 'antonini' => Antonini 9.7 tap filters
%               - 'legall'     => LeGall 5,3 tap filters.
%               - 'near_sym_a' => Near-Symmetric 5,7 tap filters.
%               - 'near_sym_b' => Near-Symmetric 13,19 tap filters.
%       par3 - qshift: 
%               - 'qshift_06' => Quarter Sample Shift Orthogonal 
%                           (Q-Shift) 10,10 tap filters, 
%                          (only 6,6 non-zero taps).
%               - 'qshift_a' =>  Q-shift 10,10 tap filters,
%                          (with 10,10 non-zero taps, unlike qshift_06).
%               - 'qshift_b' => Q-Shift 14,14 tap filters.
%               - 'qshift_c' => Q-Shift 16,16 tap filters.
%               - 'qshift_d' => Q-Shift 18,18 tap filters.
%
function X = myDTCWTSynthesis(pyr, pind, par2, par3)
       
    Ny = 4*pind(1,1);
    Nx = 4*pind(1,2);

    % Join real and imaginary parts
    M = length(pyr);
    pyr = complex(pyr(1:M/2),pyr(M/2+1:end));
    pind = pind(1:end/2,:);
    
    % Find number of levels and orientations
    norients = 6;
    nlevels = (length(pind)-1)/norients;
    
    % Reconstruct DT-CWT from pyramid
    beginIndex = 1;
    index = 1;
    for s=1:nlevels
        for w=1:norients
            endIndex = beginIndex + prod(pind(index,:)) - 1;
            Yh{s}(:,:,w) = reshape(pyr(beginIndex:endIndex),pind(index,1),pind(index,2));
            beginIndex = endIndex + 1;
            index = index + 1;
        end
    end
    Yh = Yh';
    resid = pyr(beginIndex:end);

    Nrx = Nx/2^(nlevels);
    Nry = Ny/2^(nlevels);
    
    Nrx=round(Nrx); % ---------------------------pridano ---- Michal B.
    Nry=round(Nry); % ---------------------------pridano ---- Michal B.
  
    Yl = reshape(resid,Nry,Nrx);
    
    X = dtwaveifm2(Yl,Yh,par2,par3);% It breaks if ",par5);" is included!!
return
