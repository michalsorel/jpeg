function x_est = proj2QCS(y,K,Q,size_downsampled,size_downsampled_padded,lower_bound,upper_bound,kernel,size_JPEG)            
%
%Call: x_est{m} = proj2QCS(y,K{m},Q{m},size_downsampled{m},size_downsampled_padded{m}, lower_bound{m}, upper_bound{m})
%

    right_side_resamp_pad= padarray(resample_up_down(y,size_downsampled,kernel),...
        size_downsampled_padded-size_downsampled, 'replicate', 'post'); % downsampling
    right_side_resamp_pad=Q.*bdct(right_side_resamp_pad);
    %             right_side_resamp_pad=right_side_resamp_pad - min(cat(3,max(cat(3,lower_bound{m},right_side_resamp_pad),[],3),upper_bound{m}),[],3);
    right_side_resamp_pad=right_side_resamp_pad - box_projection(right_side_resamp_pad, lower_bound, upper_bound);
    right_side_resamp_pad=K.*right_side_resamp_pad;
    right_side_resamp_pad=ibdct(right_side_resamp_pad.*Q);
    right_side_resamp_pad=right_side_resamp_pad(1:size_downsampled(1),1:size_downsampled(2)); % crop padded values
    right_side_resamp_pad=resample_up_down(right_side_resamp_pad,size_JPEG,kernel); % upsample
    x_est = y - right_side_resamp_pad;
    
    
    