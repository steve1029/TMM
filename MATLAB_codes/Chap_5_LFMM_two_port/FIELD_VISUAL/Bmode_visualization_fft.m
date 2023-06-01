% Bmode field visualization xz plane with FFT

x_inc=a/NBz;
z_inc=x_inc;
xx=[-bTx/2+x_inc/2:x_inc:bTx/2-x_inc/2];
zz=[-bTz/2+z_inc/2:z_inc:bTz/2-z_inc/2];

Nx=(length(xx)-1)/2;
Nz=(length(zz)-1)/2;

% grating region
bPG_Ey_xz=zeros(length(xx),length(zz),pcnt);
bPG_Ex_xz=zeros(length(xx),length(zz),pcnt);
bPG_Ez_xz=zeros(length(xx),length(zz),pcnt);
bPG_Hy_xz=zeros(length(xx),length(zz),pcnt);
bPG_Hx_xz=zeros(length(xx),length(zz),pcnt);
bPG_Hz_xz=zeros(length(xx),length(zz),pcnt);

bNG_Ey_xz=zeros(length(xx),length(zz),mcnt);
bNG_Ex_xz=zeros(length(xx),length(zz),mcnt);
bNG_Ez_xz=zeros(length(xx),length(zz),mcnt);
bNG_Hy_xz=zeros(length(xx),length(zz),mcnt);
bNG_Hx_xz=zeros(length(xx),length(zz),mcnt);
bNG_Hz_xz=zeros(length(xx),length(zz),mcnt);


tBpEFx=zeros(2*Nx+1,2*Nz+1,pcnt);
tBpEFy=zeros(2*Nx+1,2*Nz+1,pcnt);
tBpEFz=zeros(2*Nx+1,2*Nz+1,pcnt);
tBpHFx=zeros(2*Nx+1,2*Nz+1,pcnt);
tBpHFy=zeros(2*Nx+1,2*Nz+1,pcnt);
tBpHFz=zeros(2*Nx+1,2*Nz+1,pcnt);

tBnEFx=zeros(2*Nx+1,2*Nz+1,mcnt);
tBnEFy=zeros(2*Nx+1,2*Nz+1,mcnt);
tBnEFz=zeros(2*Nx+1,2*Nz+1,mcnt);
tBnHFx=zeros(2*Nx+1,2*Nz+1,mcnt);
tBnHFy=zeros(2*Nx+1,2*Nz+1,mcnt);
tBnHFz=zeros(2*Nx+1,2*Nz+1,mcnt);

for pbm_ind=1:pcnt
    
    for k=-nz:nz
        for m=-nx:nx
           
            tBpEFx(k+Nx+1,m+Nz+1,pbm_ind)=BpEFx(k+nz+1,1,m+nx+1,pbm_ind)*exp(j*(2*pi/bTx*k)*(bTx/2-x_inc/2-xm));
            tBpEFy(k+Nx+1,m+Nz+1,pbm_ind)=BpEFy(k+nz+1,1,m+nx+1,pbm_ind)*exp(j*(2*pi/bTx*k)*(bTx/2-x_inc/2-xm));
            tBpEFz(k+Nx+1,m+Nz+1,pbm_ind)=BpEFz(k+nz+1,1,m+nx+1,pbm_ind)*exp(j*(2*pi/bTx*k)*(bTx/2-x_inc/2-xm));
            tBpHFx(k+Nx+1,m+Nz+1,pbm_ind)=BpHFx(k+nz+1,1,m+nx+1,pbm_ind)*exp(j*(2*pi/bTx*k)*(bTx/2-x_inc/2-xm));
            tBpHFy(k+Nx+1,m+Nz+1,pbm_ind)=BpHFy(k+nz+1,1,m+nx+1,pbm_ind)*exp(j*(2*pi/bTx*k)*(bTx/2-x_inc/2-xm));
            tBpHFz(k+Nx+1,m+Nz+1,pbm_ind)=BpHFz(k+nz+1,1,m+nx+1,pbm_ind)*exp(j*(2*pi/bTx*k)*(bTx/2-x_inc/2-xm));
            
            tBnEFx(k+Nx+1,m+Nz+1,pbm_ind)=BnEFx(k+nz+1,1,m+nx+1,pbm_ind)*exp(j*(2*pi/bTx*k)*(bTx/2-x_inc/2-xm));
            tBnEFy(k+Nx+1,m+Nz+1,pbm_ind)=BnEFy(k+nz+1,1,m+nx+1,pbm_ind)*exp(j*(2*pi/bTx*k)*(bTx/2-x_inc/2-xm));
            tBnEFz(k+Nx+1,m+Nz+1,pbm_ind)=BnEFz(k+nz+1,1,m+nx+1,pbm_ind)*exp(j*(2*pi/bTx*k)*(bTx/2-x_inc/2-xm));
            tBnHFx(k+Nx+1,m+Nz+1,pbm_ind)=BnHFx(k+nz+1,1,m+nx+1,pbm_ind)*exp(j*(2*pi/bTx*k)*(bTx/2-x_inc/2-xm));
            tBnHFy(k+Nx+1,m+Nz+1,pbm_ind)=BnHFy(k+nz+1,1,m+nx+1,pbm_ind)*exp(j*(2*pi/bTx*k)*(bTx/2-x_inc/2-xm));
            tBnHFz(k+Nx+1,m+Nz+1,pbm_ind)=BnHFz(k+nz+1,1,m+nx+1,pbm_ind)*exp(j*(2*pi/bTx*k)*(bTx/2-x_inc/2-xm));
    
        end; % for m
    end; % for k
    
end; % for pbm_ind

[zG xG]=meshgrid(zz,xx); % region G
Nt=(2*Nx+1)*(2*Nz+1);

for pbm_ind=1:pcnt

    bb=odd_fftshift(tBpEFx(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    bPG_Ex_xz(:,:,pbm_ind)=bb.*exp(bp_evalue(pbm_ind)*(xG+bTx/2-x_inc/2-xm));
    
    bb=odd_fftshift(tBpEFy(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    bPG_Ey_xz(:,:,pbm_ind)=bb.*exp(bp_evalue(pbm_ind)*(xG+bTx/2-x_inc/2-xm));
    
    bb=odd_fftshift(tBpEFz(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    bPG_Ez_xz(:,:,pbm_ind)=bb.*exp(bp_evalue(pbm_ind)*(xG+bTx/2-x_inc/2-xm));
    
    bb=odd_fftshift(tBpHFx(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    bPG_Hx_xz(:,:,pbm_ind)=bb.*exp(bp_evalue(pbm_ind)*(xG+bTx/2-x_inc/2-xm));
    
    bb=odd_fftshift(tBpHFy(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    bPG_Hy_xz(:,:,pbm_ind)=bb.*exp(bp_evalue(pbm_ind)*(xG+bTx/2-x_inc/2-xm));
    
    bb=odd_fftshift(tBpHFz(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    bPG_Hz_xz(:,:,pbm_ind)=bb.*exp(bp_evalue(pbm_ind)*(xG+bTx/2-x_inc/2-xm));
    
end; % for pbm_ind=1:pcnt


for pbm_ind=1:mcnt

    bb=odd_fftshift(tBnEFx(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    bNG_Ex_xz(:,:,pbm_ind)=bb.*exp(bm_evalue(pbm_ind)*(xG+bTx/2-x_inc/2-xp));
    
    bb=odd_fftshift(tBnEFy(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    bNG_Ey_xz(:,:,pbm_ind)=bb.*exp(bm_evalue(pbm_ind)*(xG+bTx/2-x_inc/2-xp));
    
    bb=odd_fftshift(tBnEFz(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    bNG_Ez_xz(:,:,pbm_ind)=bb.*exp(bm_evalue(pbm_ind)*(xG+bTx/2-x_inc/2-xp));
    
    bb=odd_fftshift(tBnHFx(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    bNG_Hx_xz(:,:,pbm_ind)=bb.*exp(bm_evalue(pbm_ind)*(xG+bTx/2-x_inc/2-xp));
    
    bb=odd_fftshift(tBnHFy(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    bNG_Hy_xz(:,:,pbm_ind)=bb.*exp(bm_evalue(pbm_ind)*(xG+bTx/2-x_inc/2-xp));
    
    bb=odd_fftshift(tBnHFz(:,:,pbm_ind),2*Nx+1,2*Nz+1);
    cc=ifft2(bb,2*Nx+1,2*Nz+1)*Nt;
    bb=odd_ifftshift(cc,2*Nx+1,2*Nz+1);
    bNG_Hz_xz(:,:,pbm_ind)=bb.*exp(bm_evalue(pbm_ind)*(xG+bTx/2-x_inc/2-xp));
    
end; % for pbm_ind=1:pcnt








