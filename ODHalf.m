function [H im]=ODHalf(im)
%im--> input image in gray scale

im=im2double(im);
[s1 s2]=size(im);

 %Bayers 

od=[1 17 5 21 2 18 6 22;
  25 9 29 13 26 10 30 14;
  7 23 3 19 8 24 4 20;
  31 15 27 11 32 16 28 12;
  2 18 6 22 1 17 5 21;
  26 10 30 14 25 9 29 13;
  8 24 4 20 7 23 3 19;
  32 16 28 12 31 15 27 11]/32;


mask=repmat(od,round(s1/8),round(s2/8));

H=im>mask;

im=im2uint8(im);
H=im2uint8(H);

im1=imgaussfilt(im,'FilterSize',7);
H1=imgaussfilt(H,'FilterSize',7);

[peaksnr, snr] = psnr(im1, H1);
fprintf('\n The H-Peak-SNR value is %0.4f. \n', peaksnr);

% imshow(H),figure,imshow(im);

end
