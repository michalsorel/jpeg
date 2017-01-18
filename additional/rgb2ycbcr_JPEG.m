function [x, noisedev] = rgb2ycbcr_JPEG(x)
% converts RGB to YCbCr - JPEG (JFIF) standard
%
%       function [x, noisedev] = rgb2ycbcr_JPEG(x)
%
%Output:
%       x ...  image in YCbCr space
%       noisedev ... change of noise standard deviation (how many times
%                           higher than in RGB image - for each of 3 channels
%
% It expects 0-255

x = double(x);

if size(x,3)==3
    R = x(:,:,1);
    G = x(:,:,2);
    B = x(:,:,3);

    Y = 0.299 *R + 0.587 *G + 0.114 *B;
    Cb = - 0.16874* R - 0.33126 *G + 0.5 *B + 128;
    Cr = 0.5 *R - 0.41869* G - 0.08131* B + 128;

    noisedev = sqrt([0.299^2+0.587^2+0.114^2; ...
                        0.16874^2+ 0.33126^2+0.5^2; ...
                        0.5^2+0.41869^2+0.08131^2]);

    x = cat(3,Y,Cb,Cr);
end