function y=resample_up_down(x,size_y,kernel)
% simple resampling up/down for JPEG purposes

y=zeros(size_y);
switch kernel
    case 0  % pure resample
        if size(x)<size_y; % up

            y(1:2:end,1:2:end)=x;
        elseif size(x)>size_y % down
            y=x(1:2:end,1:2:end);
        else
            y=x;
        end

    case 1 % with filter
        h=[1,1,0;1,1,0;0 0 0]/4;
%         h=imrotate(h,180);        
        if size(x)>size_y; % down
            
%             y=[[y(1,1) y(1,:)];[y(:,1) y]];
            y=[[y y(:,end)];[ y(end,:) y(end,end)]];
            y=convn(x,h,'same');
            y=y(1:2:end,1:2:end);
%             y=y(1:2:end,1:2:end)/4;            
        elseif size(x)<size_y % up
            h=imrotate(h,180);
            y(1:2:end,1:2:end)=x;
%             
            y=convn(y,h,'same');
        else
            y=x;
        end
end