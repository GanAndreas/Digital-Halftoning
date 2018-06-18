function[halfpad tot_mse] = DBShalf(im)

%Applying Ordered Dithering Halftoning
im=imresize(im,[256,256]);
[s1 s2]=size(im)
H=ODHalf(im);

im=im2uint8(im);            %convert the grayscale image back to 8bit
out=im2uint8(H);          %convert the halftone image back to 8bit

%Start Dot Binary Search Algorithm
inpad=padarray(im,[1,1],'both');
halfpad1=padarray(out,[1,1],'both');
halfpad2=halfpad1;
halfpad3=halfpad1;
sum_mse(s1,s2)=0;
tot_mse=0;
t=0;
tot_mse_before=0;
dif=2;
counter=0;

while counter<=30
    while dif>=1

        for i=2:s1+1
        for j=2:s2+1
            %no mod
            block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
            temp1=halfpad2(i-1:i+1,j-1:j+1);
            block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
            temp1_gauss=imgaussfilt(temp1,1.3);
            mse1=immse(block_gauss,temp1_gauss);
            temp2=temp1;
            mse2=mse1;

            %toggle
            block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
            temp1=halfpad2(i-1:i+1,j-1:j+1);  %temp is based on the halftone image
            tog=temp1(2,2);
            if tog==255
                tog=0;
            elseif tog==0
                tog=255;
            end
            temp1(2,2)=tog;
            block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
            temp1_gauss=imgaussfilt(temp1,1.3);
            mse1=immse(block_gauss,temp1_gauss);

            if mse1<mse2
                temp2=temp1;
                mse2=mse1;
            end

            %swap 
            for y=1:3
            for x=1:3
                block=inpad(i-1:i+1,j-1:j+1);    %block is based on the grayscale image
                temp1=halfpad2(i-1:i+1,j-1:j+1);
                temp1(y,x) = halfpad2(i,j);
                temp1(2,2) = halfpad2((i-2)+y,(j-2)+x);
                block_gauss=imgaussfilt(block,1.3);    %applying gaussian filter to both block and temp
                temp1_gauss=imgaussfilt(temp1,1.3);
                mse1=immse(block_gauss,temp1_gauss);

                if mse1<mse2
                    temp2=temp1;
                    mse2=mse1;
                end
                t=t+1;
            end
            end

            sum_mse(i-1,j-1)=mse2;
            halfpad1(i-1:i+1,j-1:j+1)=temp2;

        end
        end
        tot_mse=sum(sum_mse(:))/s1^2;
        dif = abs(tot_mse_before-tot_mse);
        if counter==0
            tot_mse_before=tot_mse;
            halfpad2=halfpad1;
        else
            if tot_mse<=tot_mse_before
                tot_mse_before=tot_mse;
                halfpad2=halfpad1;
                halfpad3=halfpad1;
            else
                tot_mse_before=tot_mse_before;
                halfpad2=halfpad1;
            end
        end

            counter=counter+1;
    end
end
    
halfpad3=halfpad3(2:s1+1,2:s2+1);
imshow(im),figure,imshow(halfpad3);

im=imgaussfilt(im,1.3);
halfpad3=imgaussfilt(halfpad3,1.3);
[peaksnr, snr] = psnr(im, halfpad3);
fprintf('\n The Peak-SNR value is %0.4f. \n', peaksnr);
fprintf('\n The Total Mean Squared Error is %0.4f. \n',tot_mse);

end