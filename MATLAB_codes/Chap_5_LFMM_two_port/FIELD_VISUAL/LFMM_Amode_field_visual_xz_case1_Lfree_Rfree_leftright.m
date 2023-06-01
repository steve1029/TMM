% LFMM_Amode_field_visual_xz_case1_Lfree_Rfree_leftright

z_inc=a/NBz;
x_inc=z_inc;
zz=[0:z_inc:10*aTz-z_inc];
xx=[-aTx/2+x_inc/2:x_inc:aTx/2-x_inc/2];
zz1=zz-(10*aTz-z_inc);                      % left  freespace
zz2=zz;                                     % right freespace

% region I
Ey_xz1=zeros(length(xx),length(zz1));
Ex_xz1=zeros(length(xx),length(zz1));
Ez_xz1=zeros(length(xx),length(zz1));
Hy_xz1=zeros(length(xx),length(zz1));
Hx_xz1=zeros(length(xx),length(zz1));
Hz_xz1=zeros(length(xx),length(zz1));

% region II
Ey_xz2=zeros(length(xx),length(zz2));
Ex_xz2=zeros(length(xx),length(zz2));
Ez_xz2=zeros(length(xx),length(zz2));
Hy_xz2=zeros(length(xx),length(zz2));
Hx_xz2=zeros(length(xx),length(zz2));
Hz_xz2=zeros(length(xx),length(zz2));

y=0;
Uy=1;
Ux=0;
Einc=zeros(2*L,1);
d11=zeros(L,1);
d12=zeros(L,1);

centx=nx*NBy+ny+1;
centy=nx*NBy+ny+1;

d11(centy)=Uy;  % Kx of incident beam
d12(centx)=Ux;


Einc(1:L)=d11;
Einc(L+1:2*L)=d12;
Hinc=AVh*Einc;
Vy=Hinc(centx);
Vx=Hinc(L+centx);
Vz=-(kx0*Vx+ky0*Vy)/kz0;

ET1_temp=SAR11*Einc;
HT1_temp=-AVh*ET1_temp;
ET2_temp=SAT12*Einc;
HT2_temp=AVh*ET2_temp;

for laynt=1:SNlay

C_temp(:,laynt)=Ca(:,:,laynt)*Einc;
C_p(:,laynt)=C_temp(1:2*L,laynt);          % Positive coupling coefficients
C_n(:,laynt)=C_temp(2*L+1:4*L,laynt);      % Negative coupling coefficients

end;

ET1_y=zeros(L,1);
ET1_x=zeros(L,1);
HT1_y=zeros(L,1);
HT1_x=zeros(L,1);

ET2_y=zeros(L,1);
ET2_x=zeros(L,1);
HT2_y=zeros(L,1);
HT2_x=zeros(L,1);

ET1_y=ET1_temp(1:L);
ET1_x=ET1_temp(L+1:2*L);
HT1_y=HT1_temp(1:L);
HT1_x=HT1_temp(L+1:2*L);

ET2_y=ET2_temp(1:L);
ET2_x=ET2_temp(L+1:2*L);
HT2_y=HT2_temp(1:L);
HT2_x=HT2_temp(L+1:2*L);

ET1x_cof=zeros(NBx,NBy);
ET1y_cof=zeros(NBx,NBy);
ET1z_cof=zeros(NBx,NBy);
HT1x_cof=zeros(NBx,NBy);
HT1y_cof=zeros(NBx,NBy);
HT1z_cof=zeros(NBx,NBy);

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

      ET1y_cof(k,l)=ET1_y((k-1)*NBy+l );
      ET1x_cof(k,l)=ET1_x((k-1)*NBy+l );
	  ET1z_cof(k,l)=( Akx_vc(k)*ET1x_cof(k,l)+Aky_vc(l)*ET1y_cof(k,l) )/Akz_vc(k,l); 
      HT1y_cof(k,l)=HT1_y((k-1)*NBy+l );
      HT1x_cof(k,l)=HT1_x((k-1)*NBy+l );
	  HT1z_cof(k,l)=( Akx_vc(k)*HT1x_cof(k,l)+Aky_vc(l)*HT1y_cof(k,l) )/Akz_vc(k,l); 
   end;
end;

% diffraction efficiency

DEt1=zeros(NBx,NBy);
DEt2=zeros(NBx,NBy);
for k=1:NBx
   for l=1:NBy
      
      DEt1(k,l)=abs(abs(ET1x_cof(k,l))^2+abs(ET1y_cof(k,l))^2+abs(ET1z_cof(k,l))^2)*real(Akz_vc(k,l)/(kz0));
      DEt2(k,l)=abs(abs(ET2x_cof(k,l))^2+abs(ET2y_cof(k,l))^2+abs(ET2z_cof(k,l))^2)*real(Akz_vc(k,l)/(kz0));
      
   end;   
end;
  
total_energy=sum(sum(DEt1))+sum(sum(DEt2));

% x-z Field visualization
y=0;

% region 1
[z1 x1]=meshgrid(zz1,xx); 
  for k=1:NBx
           for l=1:NBy
           Ey_xz1=Ey_xz1+ET1y_cof(k,l)*exp(j*( Akx_vc(k)*x1+Aky_vc(l)*y-Akz_vc(k,l)*z1 ));   
           Ex_xz1=Ex_xz1+ET1x_cof(k,l)*exp(j*( Akx_vc(k)*x1+Aky_vc(l)*y-Akz_vc(k,l)*z1 ));
           Ez_xz1=Ez_xz1+ET1z_cof(k,l)*exp(j*( Akx_vc(k)*x1+Aky_vc(l)*y-Akz_vc(k,l)*z1 ));
           
           Hy_xz1=Hy_xz1+HT1y_cof(k,l)*exp(j*( Akx_vc(k)*x1+Aky_vc(l)*y-Akz_vc(k,l)*z1 ));   
           Hx_xz1=Hx_xz1+HT1x_cof(k,l)*exp(j*( Akx_vc(k)*x1+Aky_vc(l)*y-Akz_vc(k,l)*z1 ));
           Hz_xz1=Hz_xz1+HT1z_cof(k,l)*exp(j*( Akx_vc(k)*x1+Aky_vc(l)*y-Akz_vc(k,l)*z1 ));
           end;
        end;
        
       Ey_xz1=Ey_xz1+Uy*exp(j*( kx0*x1+ky0*y+kz0*z1));
       Ex_xz1=Ex_xz1+Ux*exp(j*( kx0*x1+ky0*y+kz0*z1));
       Ez_xz1=Ez_xz1+Uz*exp(j*( kx0*x1+ky0*y+kz0*z1));
       Hy_xz1=Hy_xz1+Vy*exp(j*( kx0*x1+ky0*y+kz0*z1));
       Hx_xz1=Hx_xz1+Vx*exp(j*( kx0*x1+ky0*y+kz0*z1));
       Hz_xz1=Hz_xz1+Vz*exp(j*( kx0*x1+ky0*y+kz0*z1));


% region 2

[z2 x2]=meshgrid(zz2,xx); % region G
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

zz=[0:z_inc:aTz-z_inc];
xx=[-aTx/2+x_inc/2:x_inc:aTx/2-x_inc/2];

aGEx_xz=zeros(length(xx),length(zz),SNlay);
aGEy_xz=zeros(length(xx),length(zz),SNlay);
aGEz_xz=zeros(length(xx),length(zz),SNlay);
aGHx_xz=zeros(length(xx),length(zz),SNlay);
aGHy_xz=zeros(length(xx),length(zz),SNlay);
aGHz_xz=zeros(length(xx),length(zz),SNlay);

for laynt=1:SNlay
   for pbm_ind=1:pcnt
                        
            aGEy_xz(:,:,laynt)=aGEy_xz(:,:,laynt)+C_p(pbm_ind,laynt)*aPG_Ey_xz(:,:,pbm_ind); %Sx
            aGEx_xz(:,:,laynt)=aGEx_xz(:,:,laynt)+C_p(pbm_ind,laynt)*aPG_Ex_xz(:,:,pbm_ind); %Sx
         	aGEz_xz(:,:,laynt)=aGEz_xz(:,:,laynt)+C_p(pbm_ind,laynt)*aPG_Ez_xz(:,:,pbm_ind); 
            
            aGHy_xz(:,:,laynt)=aGHy_xz(:,:,laynt)+C_p(pbm_ind,laynt)*aPG_Hy_xz(:,:,pbm_ind); %Sx
            aGHx_xz(:,:,laynt)=aGHx_xz(:,:,laynt)+C_p(pbm_ind,laynt)*aPG_Hx_xz(:,:,pbm_ind); %Sx
         	aGHz_xz(:,:,laynt)=aGHz_xz(:,:,laynt)+C_p(pbm_ind,laynt)*aPG_Hz_xz(:,:,pbm_ind); 
                  
    end;
            
      
   for pbm_ind=1:mcnt
                        
            aGEy_xz(:,:,laynt)=aGEy_xz(:,:,laynt)+C_n(pbm_ind,laynt)*aNG_Ey_xz(:,:,pbm_ind); %Sx
            aGEx_xz(:,:,laynt)=aGEx_xz(:,:,laynt)+C_n(pbm_ind,laynt)*aNG_Ex_xz(:,:,pbm_ind); %Sx
         	aGEz_xz(:,:,laynt)=aGEz_xz(:,:,laynt)+C_n(pbm_ind,laynt)*aNG_Ez_xz(:,:,pbm_ind); 
            
            aGHy_xz(:,:,laynt)=aGHy_xz(:,:,laynt)+C_n(pbm_ind,laynt)*aNG_Hy_xz(:,:,pbm_ind); %Sx
            aGHx_xz(:,:,laynt)=aGHx_xz(:,:,laynt)+C_n(pbm_ind,laynt)*aNG_Hx_xz(:,:,pbm_ind); %Sx
         	aGHz_xz(:,:,laynt)=aGHz_xz(:,:,laynt)+C_n(pbm_ind,laynt)*aNG_Hz_xz(:,:,pbm_ind); 
                  
    end;

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

figure(1); imagesc(real([Ex_xz1 AcaGEx_xz Ex_xz2])); title('Ex-xz'); colorbar;
figure(2); imagesc(real([Ey_xz1 AcaGEy_xz Ey_xz2])); title('Ey-xz'); colorbar;
figure(3); imagesc(real([Ez_xz1 AcaGEz_xz Ez_xz2])); title('Ez-xz'); colorbar;
figure(4); imagesc(real([Hx_xz1 AcaGHx_xz Hx_xz2])); title('Hx-xz'); colorbar;
figure(5); imagesc(real([Hy_xz1 AcaGHy_xz Hy_xz2])); title('Hy-xz'); colorbar;
figure(6); imagesc(real([Hz_xz1 AcaGHz_xz Hz_xz2])); title('Hz-xz'); colorbar;


for kk=0:0.2:10*2*pi
figure(7); imagesc((real([Ey_xz1 AcaGEy_xz Ey_xz2]*(exp(-i*kk)))));  caxis([-3 3]);
end;












