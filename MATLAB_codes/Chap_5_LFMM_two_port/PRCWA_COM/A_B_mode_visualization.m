% A & B mode visualization

for k=1:pcnt
mode_index=k;

figure(1); subplot(3,2,1); imagesc(real(pA_Ex_xz(:,:,mode_index))); title('pA-Ex-xz');
           subplot(3,2,2); imagesc(real(pA_Ey_xz(:,:,mode_index))); title('pA-Ey-xz');
           subplot(3,2,3); imagesc(real(pA_Ez_xz(:,:,mode_index))); title('pA-Ez-xz');
           subplot(3,2,4); imagesc(real(pA_Hx_xz(:,:,mode_index))); title('pA-Hx-xz');
           subplot(3,2,5); imagesc(real(pA_Hy_xz(:,:,mode_index))); title('pA-Hy-xz');
           subplot(3,2,6); imagesc(real(pA_Hz_xz(:,:,mode_index))); title('pA-Hz-xz');
           
           
figure(2); subplot(3,2,1); imagesc(real(pB_Ex_xz(:,:,mode_index))); title('pB-Ex-xz');
           subplot(3,2,2); imagesc(real(pB_Ey_xz(:,:,mode_index))); title('pB-Ey-xz');
           subplot(3,2,3); imagesc(real(pB_Ez_xz(:,:,mode_index))); title('pB-Ez-xz');
           subplot(3,2,4); imagesc(real(pB_Hx_xz(:,:,mode_index))); title('pB-Hx-xz');
           subplot(3,2,5); imagesc(real(pB_Hy_xz(:,:,mode_index))); title('pB-Hy-xz');
           subplot(3,2,6); imagesc(real(pB_Hz_xz(:,:,mode_index))); title('pB-Hz-xz'); 
    
           pause;
end;

for k=1:500
    
    figure(3); imagesc(real(pA_Ey_xz(:,:,16)*exp(-j*0.1*k))); title('eigenmode test');
    
end;
           