% Bmode field visualization xz plane


x_inc=a/9;
z_inc=x_inc;
xx=[0:x_inc:bTx-x_inc]';
zz=[-aTz/2+z_inc/2:z_inc:aTz/2-z_inc/2];

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


[zG xG]=meshgrid(zz,xx); % region G

y=0;

for pbm_ind=1:pcnt    
    pbm_ind
            for k=1:NBx
               for l=1:NBy
                   for m=1:NBz
              
                        
            bPG_Ex_xz(:,:,pbm_ind)=bPG_Ex_xz(:,:,pbm_ind)+BpEFx(k,l,m,pbm_ind)*exp(j*( (2*pi/bTx*(k-nx-1))*(xG-xm)+Bky_vc(l)*y+Bkz_vc(m)*zG)).*exp(bp_evalue(pbm_ind)*(xG-xm)); 
            bPG_Ey_xz(:,:,pbm_ind)=bPG_Ey_xz(:,:,pbm_ind)+BpEFy(k,l,m,pbm_ind)*exp(j*( (2*pi/bTx*(k-nx-1))*(xG-xm)+Bky_vc(l)*y+Bkz_vc(m)*zG)).*exp(bp_evalue(pbm_ind)*(xG-xm)); 
         	bPG_Ez_xz(:,:,pbm_ind)=bPG_Ez_xz(:,:,pbm_ind)+BpEFz(k,l,m,pbm_ind)*exp(j*( (2*pi/bTx*(k-nx-1))*(xG-xm)+Bky_vc(l)*y+Bkz_vc(m)*zG)).*exp(bp_evalue(pbm_ind)*(xG-xm)); 
                  
            bPG_Hx_xz(:,:,pbm_ind)=bPG_Hx_xz(:,:,pbm_ind)+BpHFx(k,l,m,pbm_ind)*exp(j*( (2*pi/bTx*(k-nx-1))*(xG-xm)+Bky_vc(l)*y+Bkz_vc(m)*zG)).*exp(bp_evalue(pbm_ind)*(xG-xm)); 
            bPG_Hy_xz(:,:,pbm_ind)=bPG_Hy_xz(:,:,pbm_ind)+BpHFy(k,l,m,pbm_ind)*exp(j*( (2*pi/bTx*(k-nx-1))*(xG-xm)+Bky_vc(l)*y+Bkz_vc(m)*zG)).*exp(bp_evalue(pbm_ind)*(xG-xm)); 
         	bPG_Hz_xz(:,:,pbm_ind)=bPG_Hz_xz(:,:,pbm_ind)+BpHFz(k,l,m,pbm_ind)*exp(j*( (2*pi/bTx*(k-nx-1))*(xG-xm)+Bky_vc(l)*y+Bkz_vc(m)*zG)).*exp(bp_evalue(pbm_ind)*(xG-xm)); 
   
                   end;
                end; % for l
      		end; % for k        
   end;

   
for pbm_ind=1:mcnt    
    pbm_ind
            for k=1:NBx
               for l=1:NBy
                   for m=1:NBz
              
                        
            bNG_Ex_xz(:,:,pbm_ind)=bNG_Ex_xz(:,:,pbm_ind)+BnEFx(k,l,m,pbm_ind)*exp(j*( (2*pi/bTx*(k-nx-1))*(xG-xm)+Bky_vc(l)*y+Bkz_vc(m)*zG)).*exp(bm_evalue(pbm_ind)*(xG-xp)); 
            bNG_Ey_xz(:,:,pbm_ind)=bNG_Ey_xz(:,:,pbm_ind)+BnEFy(k,l,m,pbm_ind)*exp(j*( (2*pi/bTx*(k-nx-1))*(xG-xm)+Bky_vc(l)*y+Bkz_vc(m)*zG)).*exp(bm_evalue(pbm_ind)*(xG-xp));
         	bNG_Ez_xz(:,:,pbm_ind)=bNG_Ez_xz(:,:,pbm_ind)+BnEFz(k,l,m,pbm_ind)*exp(j*( (2*pi/bTx*(k-nx-1))*(xG-xm)+Bky_vc(l)*y+Bkz_vc(m)*zG)).*exp(bm_evalue(pbm_ind)*(xG-xp));
   
            bNG_Hx_xz(:,:,pbm_ind)=bNG_Hx_xz(:,:,pbm_ind)+BnEFx(k,l,m,pbm_ind)*exp(j*( (2*pi/bTx*(k-nx-1))*(xG-xm)+Bky_vc(l)*y+Bkz_vc(m)*zG)).*exp(bm_evalue(pbm_ind)*(xG-xp)); 
            bNG_Hy_xz(:,:,pbm_ind)=bNG_Hy_xz(:,:,pbm_ind)+BnEFy(k,l,m,pbm_ind)*exp(j*( (2*pi/bTx*(k-nx-1))*(xG-xm)+Bky_vc(l)*y+Bkz_vc(m)*zG)).*exp(bm_evalue(pbm_ind)*(xG-xp));
         	bNG_Hz_xz(:,:,pbm_ind)=bNG_Hz_xz(:,:,pbm_ind)+BnEFz(k,l,m,pbm_ind)*exp(j*( (2*pi/bTx*(k-nx-1))*(xG-xm)+Bky_vc(l)*y+Bkz_vc(m)*zG)).*exp(bm_evalue(pbm_ind)*(xG-xp));
   
                   end;
                end; % for l
      		end; % for k        
   end;


