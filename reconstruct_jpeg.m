function [x_est,par_used] = reconstruct_jpeg(filename,varargin)
addpath('jpegtbx_1.4') % Matlab JPEG Toolbox, Phil Sallee 9/2003  <sallee@cs.ucdavis.edu>
addpath('additional')
addpath('frames')
addpath('frames/DT_CWT')

y = jpeg_read(filename);
if nargin==2;
    par_used=varargin{1};
else
    par_used=struct;
end
[x_est,snr_all,par_used] = jpeg_decoder_bregman_color_noise(y,par_used);
