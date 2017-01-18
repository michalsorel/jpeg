% example of running on jpeg image

% input
filename='images/einstein_30.jpg';
y_jpeg=imread(filename);
figure
imagesc(y_jpeg);axis image;colormap gray;colorbar vert

% reconstruct with default setting
x_est = reconstruct_jpeg(filename);
figure
imagesc(x_est);axis image;colormap gray;colorbar vert

% reconstruct with changed parameters
par.frame={'Haar'}; % available only with Wavelet Toolbox
x_est_Haar = reconstruct_jpeg(filename,par);
figure
imagesc(x_est_Haar);axis image;colormap gray;colorbar vert
