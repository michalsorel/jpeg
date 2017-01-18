function [f_thresh]=thresholding(f,q,p,threshf)
% thresholding funcion
% minimizes tau*||a||_p^p + mu/2*||f-a||^2=
%          = ||a||_p^p + q/2*||f-a||^2,      where q=mu/tau


% p = 'hard' - hard thresholding (default)
%     'soft' - soft thresholding
%     p - Lp norm

% if isempty(varargin) % third argument not passed, do default
%     type='';
%     no_of_frames=1;
% else
%     type=varargin{1};
%     no_of_frames=varargin{2};
% end

switch p
    case {0 'hard'}
        theta = (q/2)^(-1/2);
        f_thresh = hard_thresh(f,theta);
    case {1 'soft'}
        theta = 1/q; 
        f_thresh = soft_thresh(f,theta);
    otherwise % type = p, 0 < p < 1
        f_thresh = threshf(f,abs(f)); %,u_star,k)
        %warning('Thresholding not finished.');
end

end

% ---- thresholding functions -----
function f_out = hard_thresh(f,theta) % hard thresholding function    
    mask = (abs(f)>=theta);
    f_out = f.*mask;
end

function f_out = soft_thresh(f,theta) % soft thresholding function    
    f_out = sign(f).*max(abs(f)-theta,0);
end
