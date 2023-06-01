% Amode field visualization xz plane


z_inc=a/9;
x_inc=z_inc;
zz=[0:z_inc:aTz-z_inc];
xx=[-aTx/2+x_inc/2:x_inc:aTx/2-x_inc/2];

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

[zG xG]=meshgrid(zz,xx); % region G

y=0;

for pbm_ind=1:pcnt    
    pbm_ind
            for k=1:NBx
               for l=1:NBy
                   for m=1:NBz
              
                        
            aPG_Ex_xz(:,:,pbm_ind)=aPG_Ex_xz(:,:,pbm_ind)+ApEFx(k,l,m,pbm_ind)*exp(j*( Akx_vc(k)*xG+Aky_vc(l)*y+(2*pi/aTz*(m-nz-1))*(zG-zm) )).*exp(ap_evalue(pbm_ind)*(zG-zm)); %Sx
            aPG_Ey_xz(:,:,pbm_ind)=aPG_Ey_xz(:,:,pbm_ind)+ApEFy(k,l,m,pbm_ind)*exp(j*( Akx_vc(k)*xG+Aky_vc(l)*y+(2*pi/aTz*(m-nz-1))*(zG-zm) )).*exp(ap_evalue(pbm_ind)*(zG-zm)); %Sx
         	aPG_Ez_xz(:,:,pbm_ind)=aPG_Ez_xz(:,:,pbm_ind)+ApEFz(k,l,m,pbm_ind)*exp(j*( Akx_vc(k)*xG+Aky_vc(l)*y+(2*pi/aTz*(m-nz-1))*(zG-zm))).*exp(ap_evalue(pbm_ind)*(zG-zm)); 
   
            aPG_Hx_xz(:,:,pbm_ind)=aPG_Hx_xz(:,:,pbm_ind)+ApHFx(k,l,m,pbm_ind)*exp(j*( Akx_vc(k)*xG+Aky_vc(l)*y+(2*pi/aTz*(m-nz-1))*(zG-zm) )).*exp(ap_evalue(pbm_ind)*(zG-zm)); %Sx
            aPG_Hy_xz(:,:,pbm_ind)=aPG_Hy_xz(:,:,pbm_ind)+ApHFy(k,l,m,pbm_ind)*exp(j*( Akx_vc(k)*xG+Aky_vc(l)*y+(2*pi/aTz*(m-nz-1))*(zG-zm) )).*exp(ap_evalue(pbm_ind)*(zG-zm)); %Sx
         	aPG_Hz_xz(:,:,pbm_ind)=aPG_Hz_xz(:,:,pbm_ind)+ApHFz(k,l,m,pbm_ind)*exp(j*( Akx_vc(k)*xG+Aky_vc(l)*y+(2*pi/aTz*(m-nz-1))*(zG-zm))).*exp(ap_evalue(pbm_ind)*(zG-zm)); 
   
                   end;
                end; % for l
      		end; % for k        
   end;

   
for pbm_ind=1:mcnt    
    pbm_ind
            for k=1:NBx
               for l=1:NBy
                   for m=1:NBz
              
                        
            aNG_Ex_xz(:,:,pbm_ind)=aNG_Ex_xz(:,:,pbm_ind)+AnEFx(k,l,m,pbm_ind)*exp(j*( Akx_vc(k)*xG+Aky_vc(l)*y+(2*pi/aTz*(m-nz-1))*(zG-zp) )).*exp(am_evalue(pbm_ind)*(zG-zp)); %Sx
            aNG_Ey_xz(:,:,pbm_ind)=aNG_Ey_xz(:,:,pbm_ind)+AnEFy(k,l,m,pbm_ind)*exp(j*( Akx_vc(k)*xG+Aky_vc(l)*y+(2*pi/aTz*(m-nz-1))*(zG-zp) )).*exp(am_evalue(pbm_ind)*(zG-zp));
         	aNG_Ez_xz(:,:,pbm_ind)=aNG_Ez_xz(:,:,pbm_ind)+AnEFz(k,l,m,pbm_ind)*exp(j*( Akx_vc(k)*xG+Aky_vc(l)*y+(2*pi/aTz*(m-nz-1))*(zG-zp) )).*exp(am_evalue(pbm_ind)*(zG-zp));
   
            aNG_Hx_xz(:,:,pbm_ind)=aNG_Hx_xz(:,:,pbm_ind)+AnHFx(k,l,m,pbm_ind)*exp(j*( Akx_vc(k)*xG+Aky_vc(l)*y+(2*pi/aTz*(m-nz-1))*(zG-zp) )).*exp(am_evalue(pbm_ind)*(zG-zp)); %Sx
            aNG_Hy_xz(:,:,pbm_ind)=aNG_Hy_xz(:,:,pbm_ind)+AnHFy(k,l,m,pbm_ind)*exp(j*( Akx_vc(k)*xG+Aky_vc(l)*y+(2*pi/aTz*(m-nz-1))*(zG-zp) )).*exp(am_evalue(pbm_ind)*(zG-zp));
         	aNG_Hz_xz(:,:,pbm_ind)=aNG_Hz_xz(:,:,pbm_ind)+AnHFz(k,l,m,pbm_ind)*exp(j*( Akx_vc(k)*xG+Aky_vc(l)*y+(2*pi/aTz*(m-nz-1))*(zG-zp) )).*exp(am_evalue(pbm_ind)*(zG-zp));
 
                   end;
                end; % for l
      		end; % for k        
   end;
