% example of running with simulated degradation

addpath('jpegtbx_1.4\') % Matlab JPEG Toolbox, Phil Sallee 9/2003  <sallee@cs.ucdavis.edu>
addpath('images')
addpath('additional')
addpath('frames')
addpath('frames/DT_CWT')
%% Simulation
% Original image to read
soubor='hotel256'; %'hotel'; %'gantrycrane2';%'maria';%'boat'; %'einstein'
% x=imread([soubor '.jpg']);
x=imread([soubor '.png']);

% Add Gaussian noise for reconstruction with denoising
par.std_noise=0; % <0,30)
x_n=uint8(double(x)+par.std_noise*randn(size(x)));

% Set JPEG quality
par.quality=50;
soubor_jpeg=[soubor '_' num2str(par.quality) '.jpg'];
imwrite(x_n,soubor_jpeg,'Quality',par.quality)

par.original=double(x); % nove pro barvu

% Visualize
figure(1)
imagesc(x)
colormap(gray(256))
axis image
colorbar('vert')
title('Original')
caxis([0 255])
% impixelinfo
figure(2)
imagesc(x_n)
colormap(gray(256))
axis image
colorbar('vert')
title(['Original + noise: \sigma = ' num2str(par.std_noise)])
caxis([0 255])
figure(101)
x_original_ycbcr=rgb2ycbcr_JPEG(x);
for n=1:size(x_original_ycbcr,3)
    subplot(1,3,n)
    imagesc(x_original_ycbcr(:,:,n))
    colormap gray
    colorbar
    axis image
end
% impixelinfo


% Read jpeg image structure
y = jpeg_read(soubor_jpeg);
% Standard JPEG decoding
[x1 x1_ycbcr]=jpeg_decoder(y);
x1_ycbcr=double(x1_ycbcr); % to asi neni uplne korektni
% x2=imread(soubor_jpeg);
switch y.jpeg_components
    case 1
        snr_jpeg=snr(x_original_ycbcr(:,:,1),x1_ycbcr(:,:,1));
        ssim_jpeg=cal_ssim(x_original_ycbcr(:,:,1),x1_ycbcr(:,:,1),0,0);
        csnr_jpeg=csnr(x_original_ycbcr(:,:,1),x1_ycbcr(:,:,1),0,0);
    case 3
        snr_jpeg=[snr(x_original_ycbcr(:,:,1),x1_ycbcr(:,:,1)),snr(x_original_ycbcr(:,:,2),x1_ycbcr(:,:,2)),snr(x_original_ycbcr(:,:,3),x1_ycbcr(:,:,3))];
        ssim_jpeg=[cal_ssim(x_original_ycbcr(:,:,1),x1_ycbcr(:,:,1),0,0),cal_ssim(x_original_ycbcr(:,:,2),x1_ycbcr(:,:,2),0,0),cal_ssim(x_original_ycbcr(:,:,3),x1_ycbcr(:,:,3),0,0)];
        csnr_jpeg=[csnr(x_original_ycbcr(:,:,1),x1_ycbcr(:,:,1),0,0),csnr(x_original_ycbcr(:,:,2),x1_ycbcr(:,:,2),0,0),csnr(x_original_ycbcr(:,:,3),x1_ycbcr(:,:,3),0,0)];
end

% Visualize JPEG
figure(3)
imagesc(x1)
colormap(gray(256))
axis image
colorbar('vert')
title('Clasical JPEG decoder')
caxis([0 255])
% impixelinfo
figure(102)
for n=1:size(x1_ycbcr,3)
    subplot(1,3,n)
    imagesc(x1_ycbcr(:,:,n))
    colormap gray
    colorbar
    axis image
end
% impixelinfo


%% Reconstruction

% Number of iterations
par.NIter = 10; %20;%50; %
% Type of tight frame. New frame can be concetenated from basic frames.
% 'DT-CWT' - Dual-tree complex wavelets
% 'LearnedD' - Frame learned from training set
% 'Haar', 'DB2', 'DB3' - Daubechies wavelets from Wavelet Toolbox
par.frame = {'DT-CWT'}; %{'Haar'}; %{'DB2'}; %{'DB3'}; %{'LearnedD'}; %{'Haar','DB2','DB3'};
% Method used: Quantization constraint set (QCS) vs. Gaussian approximation
par.method_version = 'Gauss'; %'QCS'; %
% Joint regularization of channels
par.multichannel_regularization=true; % true - joint regularization of YCbCr channels, false - no joint regularization
% Tau is ML estimated from: JPEG or ground truth or database (tau_jpeg, tau_gt, tau_db). Default: tau_db
% par.tau_est = 'tau_db';
% Dictionary if Learned Frame
if strcmp(par.frame{1},'LearnedD')
    % load('dict16all500000.mat','dictionary'); % all Architecture database, 500000 randomly selected patches 16x16
    load('dict8all500000.mat','dictionary'); % all Architecture database, 500000 randomly selected patches 8x8
    par.dictionary=dictionary;
end
% % Learn dictionary from the original image - not available
% if strcmp(par.frame{1},'LearnedD')&&~isfield(par,'dictionary')
%     % Dictionary learning
%     % learnt_dict=learn_dictionary(double(x1_ycbcr(:,:,1)));
%     learnt_dict=learn_dictionary(double(x_original_ycbcr(:,:,1))); % z originalu
%     par.dictionary=reshape(learnt_dict,sqrt(size(learnt_dict,1)),sqrt(size(learnt_dict,1)),[]);
% end

% (Optional-advanced) Values of tau and/or mu otherwise they are set based on theory
% par.tau=0.1; % 10 - strong regularization, 0.001 - weak regularization
% par.mu= 0.1; % 0.001 - fast convergence, %0.1 - slow convergence, default: mu = tau
% par.tau= [10^-2;10^-3;10^-3];
% par.mu= [10^-2;10^-3;10^-3];

% Channel for measuring SNR
par.comp_to_measure=1;%[1 2 3];

% % Frame setup - not to be changed
% par.threshtype = {1,1,1}; % ?Type of thresholding
% pom=0;
% par.threshf = cell(numel(par.threshtype),1);
% for i = 1:numel(par.threshtype) % setup thresholding function
% %     par.threshf{i} = asetupLnormPrior(par.threshtype{i},par.tau,par.mu);
%    par.threshf{i} = asetupLnormPrior(par.threshtype{i},pom,pom);
% end

% Reconstruction
par_used=par;
[x_est,snr_all,par_used] = jpeg_decoder_bregman_color_noise(y,par_used);

% Visualize result
figure
imagesc(x_est,[0 255])
axis image
colormap(gray(256))
colorbar('vert')
title('Final estimate')
text_title=['Final estimate: ','\mu= ' num2str(par_used.mu(1)) ', \tau=' num2str(par_used.tau(1)) ...
    ', \sigma_q=' num2str(par_used.std_q(1)) ', \sigma_g=' num2str(par_used.std_noise(1))];
title(text_title)
% impixelinfo

figure
x_est_ycbcr=rgb2ycbcr_JPEG(x_est);
for n=1:y.image_components
    subplot(1,3,n)
    imagesc(x_est_ycbcr(:,:,n))
    colormap gray
    colorbar
    axis image
end
% impixelinfo

figure(106)
hold all
for n=1:length(snr_all)
    snr_s=snr_all{n};
    plot(0:(length(snr_s)-1),snr_s)
end
text_title=['\mu= ' num2str(par_used.mu(1)) ', \tau=' num2str(par_used.tau(1))];
title(text_title)
legend off
[legend_h,object_h,plot_h,text_strings] = legend('show');
text_strings((end-1):end)={ [num2str(par_used.mu(1)) ' ' num2str(par_used.tau(1)) ' ' num2str(par_used.std_q(1)) ' ' num2str(par_used.std_noise(1))] 'JPEG'};
legend(legend_h,text_strings)

% Compute improvement
switch y.jpeg_components
    case 1
        table={'','Y'};
        snr_admm=snr(x_original_ycbcr(:,:,1),x_est_ycbcr(:,:,1));
        ssim_admm=cal_ssim(x_original_ycbcr(:,:,1),x_est_ycbcr(:,:,1),0,0);
        csnr_admm=csnr(x_original_ycbcr(:,:,1),x_est_ycbcr(:,:,1),0,0);
    case 3
        table={'','Y','Cb','Cr'};
        snr_admm=[snr(x_original_ycbcr(:,:,1),x_est_ycbcr(:,:,1)),snr(x_original_ycbcr(:,:,2),x_est_ycbcr(:,:,2)),snr(x_original_ycbcr(:,:,3),x_est_ycbcr(:,:,3))];
        ssim_admm=[cal_ssim(x_original_ycbcr(:,:,1),x_est_ycbcr(:,:,1),0,0),cal_ssim(x_original_ycbcr(:,:,2),x_est_ycbcr(:,:,2),0,0),cal_ssim(x_original_ycbcr(:,:,3),x_est_ycbcr(:,:,3),0,0)];
        csnr_admm=[csnr(x_original_ycbcr(:,:,1),x_est_ycbcr(:,:,1),0,0),csnr(x_original_ycbcr(:,:,2),x_est_ycbcr(:,:,2),0,0),csnr(x_original_ycbcr(:,:,3),x_est_ycbcr(:,:,3),0,0)];
end
%% Comparison table

[table;...
    {'SNR JPEG'},num2cell(snr_jpeg);...
    {'CSNR JPEG'},num2cell(csnr_jpeg);...
    {'SSIM JPEG'},num2cell(ssim_jpeg);...
    {'ISNR ADMM'},num2cell(snr_admm-snr_jpeg);...
    {'CSNR ADMM'},num2cell(csnr_admm);...
    {'SSIM ADMM'},num2cell(ssim_admm)]