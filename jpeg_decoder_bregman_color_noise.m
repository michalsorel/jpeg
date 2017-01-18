function [x_final,snr_all par] = jpeg_decoder_bregman_color_noise(y,par)

% interp_method_up='nearest';
% interp_method_down='nearest';
fmincon_funkcional=false;
store_images=false;
kernel=1;
inimagedomain = false; %true;

groupsparsity = true;
if isfield(par,'multichannel_regularization')
    groupsparsity = par.multichannel_regularization;
end

if ~isfield(par,'comp_to_measure');   par.comp_to_measure=1;end
comp_to_measure=par.comp_to_measure;
if length(y.coef_arrays)<length(comp_to_measure)
    comp_to_measure=1:length(y.coef_arrays);
end

if ~isfield(par,'method_version'), par.method_version = 'Gaussian approximation'; end % method definition
if ~isfield(par,'frame');   par.frame={'DT-CWT'};end
if ~isfield(par,'NIter');   par.NIter=10;end

tic
% --- Initialization ---
QCSint=1;
if isfield(par,'QCS_interval'), QCSint=par.QCS_interval; end % QCS interval width

if ~isfield(par,'std_noise'), par.std_noise=0; end % noise std in original image
if length(y.coef_arrays)>1
    par.std_noise = par.std_noise.*rgb2ycbcr_noisedev; % noise decreases during rgb2ycbcr transform
else
    groupsparsity=false;
    par.multichannel_regularization=groupsparsity;
end
std_noise=par.std_noise;

[~, pom_ycbcr]=jpeg_decoder(y); % Decompress JPEG file
pom_ycbcr=double(pom_ycbcr)-128;
x_est=squeeze(mat2cell(double(pom_ycbcr),size(pom_ycbcr,1),size(pom_ycbcr,2),ones(1,size(pom_ycbcr,3))));
[analysis,synthesis,~]=frame_setup(x_est,par); % setup of analysis/synthesis functions

if isfield(par,'original')
    % gt = double(par.original)-128; % store original image RGB
    gt = double(rgb2ycbcr_JPEG(par.original))-128; % store original image YCbCr
else
    gt = pom_ycbcr;
end

if inimagedomain
    Wd = cell(1,y.image_components);
    for n=1:y.image_components
        Wd{n} = zeros(size(x_est{n}));
    end
else
    d = cell(1,y.image_components);
    for n=1:y.image_components
        %Wv{n} = doanalysis(init{n},analysis);  % Starting estimation RGB
        d{n} = doanalysis(zeros(size(x_est{n})),analysis);
    end
end

% tau estimation
if groupsparsity
    tau_jpeg = numel(pom_ycbcr(:,:,1))/jointl1norm(doanalysis(double(pom_ycbcr),analysis));
    tau_gt = numel(gt(:,:,1))/jointl1norm(doanalysis(gt,analysis));
    tau_db = 0.0546; % db Architecture DT-CWT for joint sparsity
else
    tau_jpeg = zeros(1,y.image_components);
    tau_gt = zeros(1,y.image_components);
    tau_db = [0.0572    0.4893    0.4759]; % db Architecture DT-CWT
    % tau_db=[0.0012    0.0124    0.0125]; % db Architecture LearnedD - dict8all500000.mat (not working for unknown reasons)
    for n = 1:y.image_components
        tau_jpeg(n) = numel(pom_ycbcr(:,:,1))/sum_L1_frame(double(pom_ycbcr(:,:,n)),analysis);
        tau_gt(n) = numel(gt(:,:,1))/sum_L1_frame(gt(:,:,n),analysis);
    end
end
if ~isfield(par,'tau_est'), par.tau_est='tau_db'; disp('Tau estimate based on database computation.'); end % type of tau estimator
eval(['tau_est=' par.tau_est])

disp('__________________________________________________');
disp(['MAP JPEG-denoising by ADMM, ', par.method_version])
%disp(par);
disp(['Regularization: ' par.frame{:}]);
disp(['tau_jpeg = ' num2str(tau_jpeg)]);
disp(['tau_gt = ' num2str(tau_gt)]);
disp(['tau_db = ' num2str(tau_db)]);
disp(['tau_est = ' num2str(tau_est)]);

% Prepare tau, mu, var_q
pom_ind=1:3;
pom_ind=pom_ind(isfield(par,{'tau','std_q'}));
switch sum(pom_ind)
    case 0
        disp('var_q, tau set automatically');
        var_q(1:3)=1/12;
        tau(1:3) = var_q .* tau_est;
    case 1
        disp('var_q set automatically');
        tau(1:3) = par.tau;
        var_q(1:3) = par.tau./tau_est;
    case 2
        disp('tau set automatically');
        var_q(1:3) = par.std_q.^2;
        tau(1:3) = var_q * tau_est;
    case 3
        disp(['var_q, tau taken as a parameter']);
        var_q(1:3) = par.std_q.^2;
        tau(1:3) = par.tau;
end
par.tau=tau;
par.std_q=sqrt(var_q);
if ~isfield(par,'mu'), mu_0=1*tau; else mu_0(1:3) = par.mu;disp('mu was set manually'); end % setup mu
mu= mu_0;
par.mu=mu;

% Setup thresholding
if ~isfield(par,'threshtype'),par.threshtype=repmat({1},y.image_components);end
par.threshf = cell(numel(par.threshtype),1);
for i = 1:numel(par.threshtype) % setup thresholding function
    par.threshf{i} = asetupLnormPrior(par.threshtype{i},tau,mu);
end

if store_images
    global image_iterations
    image_iterations=zeros(size(pom_ycbcr,1),size(pom_ycbcr,2),par.NIter+1,'uint8');
    image_iterations(:,:,1)=pom_ycbcr(:,:,1);
end

snr_all=cell(1,length(comp_to_measure));
sampling_factor=zeros(y.image_components,2);
k={1 1 1};

for n=1:length(comp_to_measure)
    snr_s = zeros(2,par.NIter+1);
    snr0=snr(gt(:,:,comp_to_measure(n)),double(pom_ycbcr(:,:,comp_to_measure(n)))); % YCbCr
    % snr0=snr(gt,double(pom)-128); % RGB
    snr_s(:,1) = [snr0; snr0];
    snr_s(2,:) =  snr0;
    snr_all{n}=snr_s;
end
disp_snr=cell2mat(snr_all');
initialSNRstr = ['Initial SNR = ' num2str(disp_snr(1:2:end,1)')];

size_JPEG=[y.image_height y.image_width];

for n=1:y.image_components
    % --- Precomputings for main algorithm ----
    sampling_factor(n,:)=[y.comp_info(n).v_samp_factor y.comp_info(n).h_samp_factor];
    factor=sampling_factor(n,:)./max(sampling_factor);
    if kernel
        k{n}=factor(1)*factor(2);
    end
    x_est{n} = double(pom_ycbcr(:,:,n));  % Starting estimation YCbCr
    %     x_est{n} = gt(:,:,n);  % Starting estimation YCbCr
    
    size_downsampled_padded{n}=size(y.coef_arrays{n}); % get size of downsapled channels
    size_downsampled{n}=ceil(size_JPEG.*factor);
    %     size_upsampled{n}=size_downsampled{n}./factor;
    Q{n}=repmat(1./y.quant_tables{y.comp_info(n).quant_tbl_no},size_downsampled_padded{n}/size(y.quant_tables{y.comp_info(n).quant_tbl_no},1)); % replicated quantization table over whole image
    %var_q = 1; % variance of quantization noise (hidden otherwise in tau)
    Sig=(1+k{n}*std_noise(n).^2/var_q(n)*Q{n}.^2);
    QCy=(y.coef_arrays{n}); % original DCT coefficients
    
    %Wd{n} = zeros(size(x_est{n}));
    if ~strcmp(par.method_version,'QCS') % Gaussian approximation
        CTQTQCy{n}=ibdct(Q{n}.*(1./Sig).*QCy);
        CTQTQCy{n}=CTQTQCy{n}(1:size_downsampled{n}(1),1:size_downsampled{n}(2)); % crop padded values
        CTQTQCy{n}=resample_up_down(CTQTQCy{n},size_JPEG,kernel);
        K{n}=(Q{n}.^2)./(k{n}*Q{n}.^2+mu(n).*Sig); % smart quantization table - noise version
    else
        % for QCS
        lower_bound{n} = QCy-QCSint/2;
        upper_bound{n} = QCy+QCSint/2;
        K{n}=1./(k{n}*Q{n}.^2); % smart quantization table
    end
end

% -------- funkcional used in fmincon version
if fmincon_funkcional
    f=@(x)sum(abs(x(:)));
    L1pyr=zeros(length(size_downsampled_padded),length(analysis));
    for p=1:length(size_downsampled_padded) % over image components
        for m = 1:length(analysis) % over different frames (DT-CWT,Haar)
            pyr = analysis{m}(x_est{p}(1:size_JPEG(1),1:size_JPEG(2)));
            L1pyr(p,m)=sum(cellfun(f,pyr));
        end
    end
    fmincon_funkc(1)=sum(L1pyr(:));
end
% //-------- funkcional used in fmincon version

if groupsparsity
    disp('Joint regularization (group sparsity)');
    Ax = cell(1,y.image_components);
    Axmd = cell(1,y.image_components);
else
    disp('Separate regularization for each color channel');
end

% --- Main minimization algorithm ---
for n=1:par.NIter
    x_ycbcr=[];
    if groupsparsity % color components are regularized jointly
        for m=1:y.image_components
            Ax{m} = doanalysis(x_est{m},analysis);
            Axmd{m} = subframe(Ax{m},d{m});
        end
        Wv = thresholdframejoint(Axmd,mu./tau); % BACHA na mu(m),tau(m) m je nedefinovano
    end
    for m=1:y.image_components
        if groupsparsity
            d{m} = subframe(addframe(d{m},Wv{m}), Ax{m});
        else
            if inimagedomain
                Wv{m} = threshoperator(x_est{m}-Wd{m},analysis,synthesis,mu(m)/tau(m),par);
                Wd{m} = Wd{m} - x_est{m} + Wv{m}; % tohle mozna nemusi platit poked nebude W ortogonalni
            else
                Ax = doanalysis(x_est{m},analysis);
                Wv{m} = thresholdframe(subframe(Ax,d{m}),mu(m)/tau(m),par.threshtype{m},par.threshf{m});
                d{m} = subframe(addframe(d{m},Wv{m}), Ax); % in image domain -> normal sum/subtraction
            end
        end
        
        %disp([norm(Wd(:),'fro') norm(z(:),'fro') norm(Wv(:),'fro') norm(z(:)-Wv(:),'fro') ]);
        %figure(100);imagesc(Wd,[-0.2 0.2]);pause(0.5);disp(span(Wd));
        %     if n==1
        %         figure
        %         imagesc(Wv);colormap(gray);colorbar('vert')
        %         axis image
        %         title('tresh')
        %     end
        %     snr_s(2,n+1)=snr(gt,round(Wv));
        % Step 2 in paper:
        if ~strcmp(par.method_version,'QCS')
            % a) Gaussian version
            if inimagedomain
                right_side=CTQTQCy{m}+mu(m)*(Wv{m} + Wd{m});
            else
                right_side=CTQTQCy{m}+mu(m)*dosynthesis(addframe(Wv{m},d{m}),synthesis);
            end
            
            %         CDright_side=bdct(imresize(right_side,size_downsampled{m},interp_method_down));
            right_side_resamp_pad= padarray(resample_up_down(right_side,size_downsampled{m},kernel),...
                size_downsampled_padded{m}-size_downsampled{m}, 'replicate', 'post');
            CDright_side=bdct(right_side_resamp_pad); % downsampling
            %                 figure(1000)
            %                 subplot(2,3,m)
            %                 cla
            %                 imagesc(ibdct(K{m}.*CDright_side))
            %                 colormap gray
            %                 colorbar
            %                 axis image
            
            %         DTCTdiag=imresize(ibdct(K{m}.*CDright_side),size_upsampled{n},interp_method_up);
            CTdiag=ibdct(K{m}.*CDright_side);
            CTdiag=CTdiag(1:size_downsampled{m}(1),1:size_downsampled{m}(2)); % crop padded values
            DTCTdiag=resample_up_down(CTdiag,size_JPEG,kernel); % upsample
            %         if m~=1
            %         DTCTdiag(2:2:end,:)=0;
            %         DTCTdiag(:,2:2:end)=0;
            %         end
            
            %                 subplot(2,3,3+m)
            %                 cla
            %                 imagesc(DTCTdiag)
            %                 colormap gray
            %                 colorbar
            %                 axis image
            
            x_est{m}=1/mu(m)*(right_side-DTCTdiag); % x=1/alpha*(I-DTKD)(CTQTQCy+alpha*x_est
            %        for K grayscale denoising - seems OK
            %             x_est{m}=DTCTdiag;
        else
            %     b) QCS
            %         BW-version
            %         right_side=Q.*bdct((Wv{m} + Wd{m}));
            %         %right_side=box_projection(right_side,lower_bound,upper_bound);
            %         right_side=min(cat(3,max(cat(3,lower_bound{m},right_side),[],3),upper_bound{m}),[],3);
            %         x_est{m}=ibdct(right_side./Q);
            
            % color version
            if inimagedomain
                x_est{m} = proj2QCS(Wv{m} + Wd{m},K{m},Q{m},size_downsampled{m},size_downsampled_padded{m},...
                    lower_bound{m}, upper_bound{m},kernel,size_JPEG);
            else
                x_est{m} = proj2QCS(dosynthesis(addframe(Wv{m},d{m}),synthesis),K{m},Q{m},size_downsampled{m},size_downsampled_padded{m},...
                    lower_bound{m}, upper_bound{m},kernel,size_JPEG);
            end
        end
        %                     figure(200)
        %                     imagesc(x_est{m});colormap(gray);colorbar('vert')
        %                     title(['restored, Niter=' num2str(n) ', Final SNR = ' num2str(snr(gt{m},x_est{m}))])
        
        %     x_est=round(x_est);
        %         x_est{m}(x_est{m}>127)=127;
        %         x_est{m}(x_est{m}<-128)=-128;
        x_ycbcr=cat(3,x_ycbcr,x_est{m}(1:size_JPEG(1),1:size_JPEG(2))+128);
    end
    
    %     snr_s(1,n+1)=snr(gt,round(x_est));
    %     % visualization of criterial function
    %     criterial_function_jpeg(x_est,x_est_old,y,tau,mu,par,analysis,synthesis,n)
    
    % -------- funkcional used in fmincon version
    if fmincon_funkcional
        %         f=@(x)sum(abs(x(:)));
        L1pyr=zeros(length(size_downsampled_padded),length(analysis));
        for p=1:length(size_downsampled_padded) % over image components
            for m = 1:length(analysis) % over different frames (DT-CWT,Haar)
                pyr = analysis{m}(x_est{p}(1:size_JPEG(1),1:size_JPEG(2)));
                L1pyr(p,m)=sum(cellfun(f,pyr));
            end
        end
        fmincon_funkc(n+1)=sum(L1pyr(:));
    end
    % //-------- funkcional used in fmincon version
    
    for k=1:length(comp_to_measure)
        snr_all{k}(1,n+1)=snr(gt(:,:,comp_to_measure(k)),x_ycbcr(:,:,comp_to_measure(k))-128);
    end
    
    %     figure(2000)
    %     subplot(1,2,1)
    %     cla
    %     imagesc(gt);axis image;colormap(gray(256));colorbar('vert')
    %     subplot(1,2,2)
    %     cla
    %     imagesc(x_ycbcr-128);axis image;colormap(gray(256));colorbar('vert')
    %     pause(0.1)
    
    x_final=x_ycbcr;
    if y.image_components==3
        x_final=ycbcr2rgb_JPEG(x_ycbcr);
    end
    
    %     snr_s(1,n+1)=snr(gt,x_final);
    %     x_final=uint8(x_final);
    
    
    %             figure(60)
    %             imagesc(uint8(x_final+128));axis image;
    %             colormap(gray(256))
    %             colorbar
    %             title(num2str(n))
    %                     pause
    if store_images
        image_iterations(:,:,n+1)=uint8(x_ycbcr(:,:,1));
    end
end

toc
% disp(['Final SNR = ' num2str(snr(gt,x_final))]);
% disp(['Final ISNR = ' num2str(snr(gt,x_final)-snr_s(1,1))]);

disp_snr=cell2mat(snr_all');
disp(initialSNRstr);
disp(['Final SNR = ' num2str(disp_snr(1:2:end,end)')]);
disp(['Final ISNR = ' num2str(disp_snr(1:2:end,end)'-disp_snr(1:2:end,1)')]);

x_final=uint8(x_final);

% ------ fmincon funkcional
if fmincon_funkcional
    figure(40)
    hold all
    plot(1:length(fmincon_funkc),fmincon_funkc);
    set(gca,'XScale','log')
    legend off
    [legend_h,object_h,plot_h,text_strings] = legend('show');
    text_strings(end)={ [num2str(par.mu(1)) ' ' num2str(par.tau(1)) ' ' num2str(par.std_q(1)) ' ' num2str(par.std_noise(1))]};
    %legend(legend_h,text_strings)
    legend(text_strings)
end
end
