% Amode field visualization xz plane with FFT

z_inc=a/NBz;
x_inc=z_inc;
zz=[-aTz/2+z_inc/2:z_inc:aTz/2-z_inc/2];
xx=[-aTx/2+x_inc/2:x_inc:aTx/2-x_inc/2];

Nx=(length(xx)-1)/2;
Nz=(length(zz)-1)/2;

% grating region
aPG_Ey_xz=zeros(length(xx),length(zz),pcnt);
aPG_Ex_xz=zeros(length(xx),length(zz),pcnt);
aPG_Ez_xz=zeros(length(xx),length(zz),pcnt);
aPG_Hy_xz=zeros(length(xx),length(zz),pcnt);
aPG_Hx_xz=zeros(length(xx),length(zz),pcnt);
aPG_Hz_xz=zeros(length(xx),length(zz),pcnt);

aNG_Ey_xz=zeros(length(xx),length(zz),mcnt);
aNG_Ex_xz=zeros(length(xx),length(zz),mcnt);
aNG_Ez_xz=zeros(length(xx),length(zz),mcnt);
aNG_Hy_xz=zeros(length(xx),length(zz),mcnt);
aNG_Hx_xz=zeros(length(xx),length(zz),mcnt);
aNG_Hz_xz=zeros(length(xx),length(zz),mcnt);

tApEFx=zeros(2*Nx+1,2*Nz+1,pcnt);
tApEFy=zeros(2*Nx+1,2*Nz+1,pcnt);
tApEFz=zeros(2*Nx+1,2*Nz+1,pcnt);
tApHFx=zeros(2*Nx+1,2*Nz+1,pcnt);
tApHFy=zeros(2*Nx+1,2*Nz+1,pcnt);
tApHFz=zeros(2*Nx+1,2*Nz+1,pcnt);

tAnEFx=zeros(2*Nx+1,2*Nz+1,mcnt);
tAnEFy=zeros(2*Nx+1,2*Nz+1,mcnt);
tAnEFz=zeros(2*Nx+1,2*Nz+1,mcnt);
tAnHFx=zeros(2*Nx+1,2*Nz+1,mcnt);
tAnHFy=zeros(2*Nx+1,2*Nz+1,mcnt);
tAnHFz=zeros(2*Nx+1,2*Nz+1,mcnt);


for pbm_ind=1:pcnt
    
        for k=-nx:nx
            for m=-nz:nz

            tApEFx(k+Nx+1,m+Nz+1,pbm_ind)=ApEFx(k+nx+1,1,m+nz+1,pbm_ind)*exp(j*(2*pi/aTz*m)*(aTz/2-z_inc/2-zm));
            tApEFy(k+Nx+1,m+Nz+1,pbm_ind)=ApEFy(k+nx+1,1,m+nz+1,pbm_ind)*exp(j*(2*pi/aTz*m)*(aTz/2-z_inc/2-zm));
            tApEFz(k+Nx+1,m+Nz+1,pbm_ind)=ApEFz(k+nx+1,1,m+nz+1,pbm_ind)*exp(j*(2*pi/aTz*m)*(aTz/2-z_inc/2-zm));
            tApHFx(k+Nx+1,m+Nz+1,pbm_ind)=ApHFx(k+nx+1,1,m+nz+1,pbm_ind)*exp(j*(2*pi/aTz*m)*(aTz/2-z_inc/2-zm));
            tApHFy(k+Nx+1,m+Nz+1,pbm_ind)=ApHFy(k+nx+1,1,m+nz+1,pbm_ind)*exp(j*(2*pi/aTz*m)*(aTz/2-z_inc/2-zm));
            tApHFz(k+Nx+1,m+Nz+1,pbm_ind)=ApHFz(k+nx+1,1,m+nz+1,pbm_ind)*exp(j*(2*pi/aTz*m)*(aTz/2-z_inc/2-zm));
            
            tAnEFx(k+Nx+1,m+Nz+1,pbm_ind)=AnEFx(k+nx+1,1,m+nz+1,pbm_ind)*exp(j*(2*pi/aTz*m)*(aTz/2-z_inc/2-zm));
            tAnEFy(k+Nx+1,m+Nz+1,pbm_ind)=AnEFy(k+nx+1,1,m+nz+1,pbm_ind)*exp(j*(2*pi/aTz*m)*(aTz/2-z_inc/2-zm));
            tAnEFz(k+Nx+1,m+Nz+1,pbm_ind)=AnEFz(k+nx+1,1,m+nz+1,pbm_ind)*exp(j*(2*pi/aTz*m)*(aTz/2-z_inc/2-zm));
            tAnHFx(k+Nx+1,m+Nz+1,pbm_ind)=AnHFx(k+nx+1,1,m+nz+1,pbm_ind)*exp(j*(2*pi/aTz*m)*(aTz/2-z_inc/2-zm));
            tAnHFy(k+Nx+1,m+Nz+1,pbm_ind)=AnHFy(k+nx+1,1,m+nz+1,pbm_ind)*exp(j*(2*pi/aTz*m)*(aTz/2-z_inc/2-zm));
            tAnHFz(k+Nx+1,m+Nz+1,pbm_ind)=AnHFz(k+nx+1,1,m+nz+1,pbm_ind)*exp(j*(2*pi/aTz*m)*(aTz/2-z_inc/2-zm));
            
            end;    % for m
        end;    % for k

end; % for pbm_ind


[zG xG]=meshgrid(zz,xx); % region G
Nt=(2*Nx+1)*(2*Nz+1);
for pbm_ind=1:pcnt

    bb=odd_fftshift(tApEFx(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    aPG_Ex_xz(:,:,pbm_ind)=bb.*exp(ap_evalue(pbm_ind)*(zG+aTz/2-z_inc/2-zm));
    
    bb=odd_fftshift(tApEFy(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    aPG_Ey_xz(:,:,pbm_ind)=bb.*exp(ap_evalue(pbm_ind)*(zG+aTz/2-z_inc/2-zm));
    
    bb=odd_fftshift(tApEFz(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    aPG_Ez_xz(:,:,pbm_ind)=bb.*exp(ap_evalue(pbm_ind)*(zG+aTz/2-z_inc/2-zm));
    
    bb=odd_fftshift(tApHFx(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    aPG_Hx_xz(:,:,pbm_ind)=bb.*exp(ap_evalue(pbm_ind)*(zG+aTz/2-z_inc/2-zm));
    
    bb=odd_fftshift(tApHFy(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    aPG_Hy_xz(:,:,pbm_ind)=bb.*exp(ap_evalue(pbm_ind)*(zG+aTz/2-z_inc/2-zm));
    
    bb=odd_fftshift(tApHFz(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    aPG_Hz_xz(:,:,pbm_ind)=bb.*exp(ap_evalue(pbm_ind)*(zG+aTz/2-z_inc/2-zm));
    
end; % for pbm_ind=1:pcnt


for pbm_ind=1:mcnt
    
    bb=odd_fftshift(tAnEFx(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    aNG_Ex_xz(:,:,pbm_ind)=bb.*exp(am_evalue(pbm_ind)*(zG+aTz/2-z_inc/2-zp));
    
    bb=odd_fftshift(tAnEFy(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    aNG_Ey_xz(:,:,pbm_ind)=bb.*exp(am_evalue(pbm_ind)*(zG+aTz/2-z_inc/2-zp));
    
    bb=odd_fftshift(tAnEFz(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    aNG_Ez_xz(:,:,pbm_ind)=bb.*exp(am_evalue(pbm_ind)*(zG+aTz/2-z_inc/2-zp));
    
    bb=odd_fftshift(tAnHFx(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    aNG_Hx_xz(:,:,pbm_ind)=bb.*exp(am_evalue(pbm_ind)*(zG+aTz/2-z_inc/2-zp));
    
    bb=odd_fftshift(tAnHFy(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    aNG_Hy_xz(:,:,pbm_ind)=bb.*exp(am_evalue(pbm_ind)*(zG+aTz/2-z_inc/2-zp));
    
    bb=odd_fftshift(tAnHFz(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    aNG_Hz_xz(:,:,pbm_ind)=bb.*exp(am_evalue(pbm_ind)*(zG+aTz/2-z_inc/2-zp));
    
    
end; % for pbm_ind=1:mcnt

            
            
