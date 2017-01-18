function noisedev = rgb2ycbcr_noisedev
%rgb2ycbcr_noisedev  attenuation of noise standard deviation during conversion RGB to YCbCr 
%
%       function noisedev = rgb2ycbcr_noisedev
%
%Output:
%       noisedev ... change of noise standard deviation (how many times
%                           higher than in RGB image - for each of 3 channels
%
%See also rgb2ycbcr and rgb2ycbcr_JPEG
%
%Michal Sorel 2016

%Y = 0.299 *R + 0.587 *G + 0.114 *B;
%Cb = - 0.16874* R - 0.33126 *G + 0.5 *B + 128;
%Cr = 0.5 *R - 0.41869* G - 0.08131* B + 128;

noisedev = sqrt([0.299^2+0.587^2+0.114^2; ...
                    0.16874^2+ 0.33126^2+0.5^2; ...
                    0.5^2+0.41869^2+0.08131^2]);

