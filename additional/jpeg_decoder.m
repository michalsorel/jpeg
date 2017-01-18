function [x x_ycbcr]=jpeg_decoder(y)
% simple jpeg decoder //
% x=jpeg_decoder(jpeg_object)
kernel=1;
for n=1:y.image_components
    x_help{n}=ibdct(dequantize(y.coef_arrays{n},y.quant_tables{y.comp_info(n).quant_tbl_no}));
sampling_factor(n,:)=[y.comp_info(n).v_samp_factor y.comp_info(n).h_samp_factor];
end

x_ycbcr=zeros([y.image_height,y.image_width,y.image_components]);
for n=1:y.image_components
    factor=sampling_factor(n,:)./max(sampling_factor);
%         x_ycbcr(:,:,n) = imresize(x_help{n},[y.image_height y.image_width]./factor,'bilinear');
%         x_help{n} = imresize(x_help{n},size(x_help{n})./factor,'bilinear');
        x_help{n} = resample_up_down(x_help{n},size(x_help{n})./factor,kernel)/prod(factor);

        x_ycbcr(:,:,n) = x_help{n}(1:y.image_height,1:y.image_width); 
%     [rows cols]=size(x_help{n});
%     step=[y.image_height/rows,y.image_width/cols];
%     off=step/2-0.5;
%     temp=x_help{n};
%     temp=[temp(1,1) temp(1,:) temp(1,end); temp(:,1) temp temp(:,end); temp(end,1) temp(end,:) temp(end,end)];
%     [X,Y] = meshgrid((1-step(1)+off(1)):step(1):(y.image_height+step(1)),(1-step(2)+off(2)):step(2):(y.image_width+step(2)));
%     [XI,YI] = meshgrid(1:y.image_height,1:y.image_width);
%     x_ycbcr(:,:,n)  = interp2(X,Y,temp,XI,YI);

%         figure(4)
%         subplot(2,3,n)
%         imagesc(x_ycbcr(:,:,n))
%         colormap gray
%         colorbar
%         axis image
end

% x=int8(x);
% x=uint8(x);
% x = im2uint8((x+128)/255); % 128 kvuli offsetu pro uint8
% x = uint8((x+128)); % 128 kvuli offsetu pro uint8
% x = ((x+128)/255); % 128 kvuli offsetu pro uint8
% x=int8(x);

x=double(uint8(x_ycbcr+128));


if y.image_components==3
    %     x=ycbcr2rgb(x);
x=ycbcr2rgb_JPEG(x);
%     Y=x(:,:,1);
%     cb=x(:,:,2);
%     cr=x(:,:,3);
%     r=Y                 +1.402  *(cr-128);
%     g=Y-0.34414*(cb-128)-0.71414*(cr-128);
%     b=Y+1.772  *(cb-128);
%     x=cat(3,r,g,b);
end
x=uint8(x);
x_ycbcr=uint8(x_ycbcr+128);
% x = im2uint8(x);
end