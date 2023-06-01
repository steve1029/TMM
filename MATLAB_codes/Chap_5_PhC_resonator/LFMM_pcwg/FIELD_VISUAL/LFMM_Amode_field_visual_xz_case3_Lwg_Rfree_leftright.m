% LFMM_Amode_field_visual_xz_case3_Lwg_Rfree_leftright

D=zeros(2*L,1);
inc_mode=33;
D(inc_mode)=1;

C_n=ALR11*D;         % negative coupling coefficient
ET2_temp=ALT12*D;    % positive coupling coefficient
HT2_temp=AVh*ET2_temp;

ET2_y=zeros(L,1);
ET2_x=zeros(L,1);
HT2_y=zeros(L,1);
HT2_x=zeros(L,1);

ET2_y=ET2_temp(1:L);
ET2_x=ET2_temp(L+1:2*L);
HT2_y=HT2_temp(1:L);
HT2_x=HT2_temp(L+1:2*L);

ET2x_cof=zeros(NBx,NBy);
ET2y_cof=zeros(NBx,NBy);
ET2z_cof=zeros(NBx,NBy);
HT2x_cof=zeros(NBx,NBy);
HT2y_cof=zeros(NBx,NBy);
HT2z_cof=zeros(NBx,NBy);

for k=1:NBx
   for l=1:NBy
      
      ET2y_cof(k,l)=ET2_y((k-1)*NBy+l);
      ET2x_cof(k,l)=ET2_x((k-1)*NBy+l );
      ET2z_cof(k,l)=-( Akx_vc(k)*ET2x_cof(k,l)+Aky_vc(l)*ET2y_cof(k,l) )/Akz_vc(k,l); 
      HT2y_cof(k,l)=HT2_y((k-1)*NBy+l);
      HT2x_cof(k,l)=HT2_x((k-1)*NBy+l );
      HT2z_cof(k,l)=-( Akx_vc(k)*HT2x_cof(k,l)+Aky_vc(l)*HT2y_cof(k,l) )/Akz_vc(k,l); 
 
   end;
end;


SNlay=60;
% x-z Field visualization
y=0;
z_inc=a/NBz;
x_inc=z_inc;
zz=[0:z_inc:SNlay*aTz-z_inc];
xx=[-aTx/2+x_inc/2:x_inc:aTx/2-x_inc/2];
zz1=zz-(SNlay*aTz-z_inc);                      % left  space
zz2=zz;                                        % right space

% region 2
Ey_xz2=zeros(length(xx),length(zz2));
Ex_xz2=zeros(length(xx),length(zz2));
Ez_xz2=zeros(length(xx),length(zz2));
Hy_xz2=zeros(length(xx),length(zz2));
Hx_xz2=zeros(length(xx),length(zz2));
Hz_xz2=zeros(length(xx),length(zz2));


% region 2
[z2 x2]=meshgrid(zz2,xx); 
 y=0;
    for k=1:NBx
         for l=1:NBy
           Ey_xz2=Ey_xz2+ET2y_cof(k,l)*exp(j*Akz_vc(k,l)*z2).*exp(j*( Akx_vc(k)*x2+Aky_vc(l)*y ));   
           Ex_xz2=Ex_xz2+ET2x_cof(k,l)*exp(j*Akz_vc(k,l)*z2).*exp(j*( Akx_vc(k)*x2+Aky_vc(l)*y ));
           Ez_xz2=Ez_xz2+ET2z_cof(k,l)*exp(j*Akz_vc(k,l)*z2).*exp(j*( Akx_vc(k)*x2+Aky_vc(l)*y ));
           
           Hy_xz2=Hy_xz2+HT2y_cof(k,l)*exp(j*Akz_vc(k,l)*z2).*exp(j*( Akx_vc(k)*x2+Aky_vc(l)*y ));   
           Hx_xz2=Hx_xz2+HT2x_cof(k,l)*exp(j*Akz_vc(k,l)*z2).*exp(j*( Akx_vc(k)*x2+Aky_vc(l)*y ));
           Hz_xz2=Hz_xz2+HT2z_cof(k,l)*exp(j*Akz_vc(k,l)*z2).*exp(j*( Akx_vc(k)*x2+Aky_vc(l)*y ));
           end;
        end;
 
      
% Grating region

[xlen zlen]=size(aNG_Ey_xz(:,:,1));  
aGEx_xz=zeros(xlen,zlen,SNlay);
aGEy_xz=zeros(xlen,zlen,SNlay);
aGEz_xz=zeros(xlen,zlen,SNlay);
aGHx_xz=zeros(xlen,zlen,SNlay);
aGHy_xz=zeros(xlen,zlen,SNlay);
aGHz_xz=zeros(xlen,zlen,SNlay);

for laynt=1:SNlay
   for pbm_ind=1:mcnt
                        
            aGEy_xz(:,:,laynt)=aGEy_xz(:,:,laynt)+C_n(pbm_ind)*aNG_Ey_xz(:,:,pbm_ind)*exp(-am_evalue(pbm_ind)*aTz*(SNlay-laynt)); %Sx
            aGEx_xz(:,:,laynt)=aGEx_xz(:,:,laynt)+C_n(pbm_ind)*aNG_Ex_xz(:,:,pbm_ind)*exp(-am_evalue(pbm_ind)*aTz*(SNlay-laynt)); %Sx
         	aGEz_xz(:,:,laynt)=aGEz_xz(:,:,laynt)+C_n(pbm_ind)*aNG_Ez_xz(:,:,pbm_ind)*exp(-am_evalue(pbm_ind)*aTz*(SNlay-laynt)); 
            
            aGHy_xz(:,:,laynt)=aGHy_xz(:,:,laynt)+C_n(pbm_ind)*aNG_Hy_xz(:,:,pbm_ind)*exp(-am_evalue(pbm_ind)*aTz*(SNlay-laynt)); %Sx
            aGHx_xz(:,:,laynt)=aGHx_xz(:,:,laynt)+C_n(pbm_ind)*aNG_Hx_xz(:,:,pbm_ind)*exp(-am_evalue(pbm_ind)*aTz*(SNlay-laynt)); %Sx
         	aGHz_xz(:,:,laynt)=aGHz_xz(:,:,laynt)+C_n(pbm_ind)*aNG_Hz_xz(:,:,pbm_ind)*exp(-am_evalue(pbm_ind)*aTz*(SNlay-laynt)); 
                  
    end;
            % incidence mode
            aGEy_xz(:,:,laynt)=aGEy_xz(:,:,laynt)+aPG_Ey_xz(:,:,inc_mode)*exp(-ap_evalue(inc_mode)*aTz*(SNlay-laynt+1)); %Sx
            aGEx_xz(:,:,laynt)=aGEx_xz(:,:,laynt)+aPG_Ex_xz(:,:,inc_mode)*exp(-ap_evalue(inc_mode)*aTz*(SNlay-laynt+1)); %Sx
         	aGEz_xz(:,:,laynt)=aGEz_xz(:,:,laynt)+aPG_Ez_xz(:,:,inc_mode)*exp(-ap_evalue(inc_mode)*aTz*(SNlay-laynt+1)); 
            
            aGHy_xz(:,:,laynt)=aGHy_xz(:,:,laynt)+aPG_Hy_xz(:,:,inc_mode)*exp(-ap_evalue(inc_mode)*aTz*(SNlay-laynt+1)); %Sx
            aGHx_xz(:,:,laynt)=aGHx_xz(:,:,laynt)+aPG_Hx_xz(:,:,inc_mode)*exp(-ap_evalue(inc_mode)*aTz*(SNlay-laynt+1)); %Sx
         	aGHz_xz(:,:,laynt)=aGHz_xz(:,:,laynt)+aPG_Hz_xz(:,:,inc_mode)*exp(-ap_evalue(inc_mode)*aTz*(SNlay-laynt+1)); 
           
    
end; % for laynt


AcaGEx_xz=aGEx_xz(:,:,1);
AcaGEy_xz=aGEy_xz(:,:,1);
AcaGEz_xz=aGEz_xz(:,:,1);

AcaGHx_xz=aGHx_xz(:,:,1);
AcaGHy_xz=aGHy_xz(:,:,1);
AcaGHz_xz=aGHz_xz(:,:,1);

for laynt=2:SNlay
    AcaGEx_xz=[AcaGEx_xz aGEx_xz(:,:,laynt)];
    AcaGEy_xz=[AcaGEy_xz aGEy_xz(:,:,laynt)];
    AcaGEz_xz=[AcaGEz_xz aGEz_xz(:,:,laynt)];
    
    AcaGHx_xz=[AcaGHx_xz aGHx_xz(:,:,laynt)];
    AcaGHy_xz=[AcaGHy_xz aGHy_xz(:,:,laynt)];
    AcaGHz_xz=[AcaGHz_xz aGHz_xz(:,:,laynt)];
end;


figure(1); imagesc(real([AcaGEx_xz Ex_xz2])); colorbar;  title('Ex-xz');
figure(2); imagesc(real([AcaGEy_xz Ey_xz2])); colorbar;  title('Ey-xz');
figure(3); imagesc(real([AcaGEz_xz Ez_xz2])); colorbar;  title('Ez-xz'); 

figure(4); imagesc(real([AcaGHx_xz Hx_xz2])); colorbar;  title('Hx-xz'); 
figure(5); imagesc(real([AcaGHy_xz Hy_xz2])); colorbar;  title('Hy-xz');
figure(6); imagesc(real([AcaGHz_xz Hz_xz2])); colorbar;  title('Hz-xz');











