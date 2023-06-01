function bb=DFT_conv(gy,NBz)


        %bb=odd_fftshift(gy,1,NBz);
        bb=gy;
        cc=fft2(bb,1,NBz);
        bb=odd_ifftshift(cc,1,NBz);
        bb=bb/NBz;
        
%             cc=odd_fftshift(bb,1,NBz);
%             bb=ifft2(cc,1,NBz);
%             bb=odd_ifftshift(bb,1,NBz)*NBz;
%             
        