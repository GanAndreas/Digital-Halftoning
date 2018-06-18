function[halfpad,tot_mse] = DBShalf(im)

%Applying Ordered Dithering Halftoning
im=rgb2gray(im);
im=imresize(im,[256,256]);
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
out=im>mask;

im=im2uint8(im);            %convert the grayscale image back to 8bit
out=im2uint8(out);          %convert the halftone image back to 8bit

%Start Dot Binary Search Algorithm
inpad=padarray(im,[1,1],'both');
halfpad=padarray(out,[1,1],'both');
mse(10)=0;
sum_mse(s1,s2)=0;
tot_mse=0;
tot_mse_before=0;
dif=2;
counter=0;

while dif>=0.1
    
    for i=2:s1+1
    for j=2:s2+1
        %no mod
        block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
        temp1=halfpad(i-1:i+1,j-1:j+1);
        tot_sum = 0;
        block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
        temp1_gauss=imgaussfilt(temp1,1.3);
%         for y=1:3
%         for x=1:3
%             sum_1=abs(block_gauss(y,x)-temp1_gauss(y,x));
%             tot_sum = double(tot_sum + sum_1);
%         end
%         end
        err=immse(block_gauss,temp1_gauss);
        mse(1)=err;

        %toggle
        block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
        temp2=halfpad(i-1:i+1,j-1:j+1);  %temp is based on the halftone image
        tog=temp2(2,2);
        if tog==255
            tog=0;
        elseif tog==0
            tog=255;
        end
        temp2(2,2)=tog;
        tot_sum = 0;
        block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
        temp2_gauss=imgaussfilt(temp2,1.3);
%         for y=1:3
%         for x=1:3
%             sum_1=abs(block_gauss(y,x)-temp2_gauss(y,x));
%             tot_sum = double(tot_sum + sum_1);
%         end
%         end
        err=immse(block_gauss,temp2_gauss);
        mse(2)=err;

        %swap 1
        block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
        temp3=halfpad(i-1:i+1,j-1:j+1);  %temp is based on the halftone image
        temp3(1,1) = halfpad(i,j);
        temp3(2,2) = halfpad(i-1,j-1);
        tot_sum = 0;
        block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
        temp3_gauss=imgaussfilt(temp3,1.3);
%         for y=1:3
%         for x=1:3
%             sum_1=abs(block_gauss(y,x)-temp3_gauss(y,x));
%             tot_sum = double(tot_sum + sum_1);
%         end
%         end
        err=immse(block_gauss,temp3_gauss);
        mse(3)=err;

        %swap 2
        block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
        temp4=halfpad(i-1:i+1,j-1:j+1);
        temp4(1,2) = halfpad(i,j);
        temp4(2,2) = halfpad(i-1,j);
        tot_sum = 0;
        block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
        temp4_gauss=imgaussfilt(temp4,1.3);
        for y=1:3
        for x=1:3
            sum_1=abs(block_gauss(y,x)-temp4_gauss(y,x));
            tot_sum = double(tot_sum + sum_1);
        end
        end
        mse(4)=tot_sum/9;

        %swap 3
        block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
        temp5=halfpad(i-1:i+1,j-1:j+1);
        temp5(1,3) = halfpad(i,j);
        temp5(2,2) = halfpad(i-1,j+1);
        tot_sum = 0;
        block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
        temp5_gauss=imgaussfilt(temp5,1.3);
%         for y=1:3
%         for x=1:3
%             sum_1=abs(block_gauss(y,x)-temp5_gauss(y,x));
%             tot_sum = double(tot_sum + sum_1);
%         end
%         end
        err=immse(block_gauss,temp5_gauss);
        mse(5)=err;

        %swap 4
        block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
        temp6=halfpad(i-1:i+1,j-1:j+1);
        temp6(2,1) = halfpad(i,j);
        temp6(2,2) = halfpad(i,j-1);
        tot_sum = 0;
        block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
        temp6_gauss=imgaussfilt(temp6,1.3);
%         for y=1:3
%         for x=1:3
%             sum_1=abs(block_gauss(y,x)-temp6_gauss(y,x));
%             tot_sum = double(tot_sum + sum_1);
%         end
%         end
        err=immse(block_gauss,temp6_gauss);
        mse(6)=err;

        %swap 5
        block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
        temp7=halfpad(i-1:i+1,j-1:j+1);
        temp7(2,3) = halfpad(i,j);
        temp7(2,2) = halfpad(i,j+1);
        tot_sum = 0;
        block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
        temp7_gauss=imgaussfilt(temp7,1.3);
%         for y=1:3
%         for x=1:3
%             sum_1=abs(block_gauss(y,x)-temp7_gauss(y,x));
%             tot_sum = double(tot_sum + sum_1);
%         end
%         end
        err=immse(block_gauss,temp7_gauss);
        mse(7)=err;

        %swap 6
        block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
        temp8=halfpad(i-1:i+1,j-1:j+1);
        temp8(3,1) = halfpad(i,j);
        temp8(2,2) = halfpad(i+1,j-1);
        tot_sum = 0;
        block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
        temp8_gauss=imgaussfilt(temp8,1.3);
%         for y=1:3
%         for x=1:3
%             sum_1=abs(block_gauss(y,x)-temp8_gauss(y,x));
%             tot_sum = double(tot_sum + sum_1);
%         end
%         end
        err=immse(block_gauss,temp8_gauss);
        mse(8)=err;

        %swap 7
        block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
        temp9=halfpad(i-1:i+1,j-1:j+1);
        temp9(3,2) = halfpad(i,j);
        temp9(2,2) = halfpad(i+1,j);
        tot_sum = 0;
        block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
        temp9_gauss=imgaussfilt(temp9,1.3);
%         for y=1:3
%         for x=1:3
%             sum_1=abs(block_gauss(y,x)-temp9_gauss(y,x));
%             tot_sum = double(tot_sum + sum_1);
%         end
%         end
        err=immse(block_gauss,temp9_gauss);
        mse(9)=err;

        %swap 8
        block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
        temp10=halfpad(i-1:i+1,j-1:j+1);
        temp10(2,3) = halfpad(i,j);
        temp10(2,2) = halfpad(i+1,j+1);
        tot_sum = 0;
        block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
        temp10_gauss=imgaussfilt(temp10,1.3);
%         for y=1:3
%         for x=1:3
%             sum_1=abs(block_gauss(y,x)-temp9_gauss(y,x));
%             tot_sum = double(tot_sum + sum_1);
%         end
%         end
        err=immse(block_gauss,temp10_gauss);
        mse(10)=err;

        mse1=min(mse(:));               %find the smallest MSE
        p = find(mse==mse1);            
        p = min(p(:));
        sum_mse(i-1,j-1)=mse1;

        switch p
            case 1
                halfpad(i-1:i+1,j-1:j+1)=temp1;
            case 2
                halfpad(i-1:i+1,j-1:j+1)=temp2;
            case 3
                halfpad(i-1:i+1,j-1:j+1)=temp3;
            case 4
                halfpad(i-1:i+1,j-1:j+1)=temp4;
            case 5
                halfpad(i-1:i+1,j-1:j+1)=temp5;
            case 6
                halfpad(i-1:i+1,j-1:j+1)=temp6;
            case 7
                halfpad(i-1:i+1,j-1:j+1)=temp7;
            case 8
                halfpad(i-1:i+1,j-1:j+1)=temp7;
            case 9
                halfpad(i-1:i+1,j-1:j+1)=temp9;
            case 10
                halfpad(i-1:i+1,j-1:j+1)=temp10;
        end

    end    
    end
    tot_mse=sum(sum_mse(:))/s1^2;
    dif = abs(tot_mse_before-tot_mse);
    if counter==0
        tot_mse_before=tot_mse;
%         dif=1;
    else
        if tot_mse<=tot_mse_before
            tot_mse_before=tot_mse;
        else
            tot_mse_before=tot_mse_before;
        end
    end

        counter=counter+1;
end    

halfpad=halfpad(2:s1+1,2:s2+1);
imshow(im),figure,imshow(halfpad);

im=imgaussfilt(im,1.3);
halfpad=imgaussfilt(halfpad,1.3);
[peaksnr, snr] = psnr(im, halfpad);
fprintf('\n The Peak-SNR value is %0.4f. \n', peaksnr);
fprintf('\n The Total Mean Squared Error is %0.4f. \n',tot_mse);

end