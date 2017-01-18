function x=ycbcr2rgb_JPEG(x)
% converts YCbCr to RGB - JPEG (JFIF) standard
x = double(x);
Y=x(:,:,1);
cb=x(:,:,2);
cr=x(:,:,3);

r=Y                 +1.402  *(cr-128);
g=Y-0.34414*(cb-128)-0.71414*(cr-128);
b=Y+1.772  *(cb-128);

x=cat(3,r,g,b);