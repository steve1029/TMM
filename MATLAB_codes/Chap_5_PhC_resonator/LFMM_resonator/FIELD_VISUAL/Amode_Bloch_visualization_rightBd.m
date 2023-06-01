% Amode field visualization in x-z plane
% Right grating
% freespace - grating

Uy=1
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

ET1_temp=ARR11*Einc;
HT1_temp=-AVh*ET1_temp;

C_p=ART12*Einc; % positive coupling coefficient

ET1_y=zeros(L,1);
ET1_x=zeros(L,1);
HT1_y=zeros(L,1);
HT1_x=zeros(L,1);

ET1_y=ET1_temp(1:L);
ET1_x=ET1_temp(L+1:2*L);
HT1_y=HT1_temp(1:L);
HT1_x=HT1_temp(L+1:2*L);

ET1x_cof=zeros(NBx,NBy);
ET1y_cof=zeros(NBx,NBy);
ET1z_cof=zeros(NBx,NBy);
HT1x_cof=zeros(NBx,NBy);
HT1y_cof=zeros(NBx,NBy);
HT1z_cof=zeros(NBx,NBy);


for k=1:NBx
   for l=1:NBy
      % electric field
      ET1y_cof(k,l)=ET1_y((k-1)*NBy+l );
      ET1x_cof(k,l)=ET1_x((k-1)*NBy+l );
	  ET1z_cof(k,l)=( Akx_vc(k)*ET1x_cof(k,l)+Aky_vc(l)*ET1y_cof(k,l) )/Akz_vc(k,l); 
      % magnetic field
      HT1y_cof(k,l)=HT1_y((k-1)*NBy+l );
      HT1x_cof(k,l)=HT1_x((k-1)*NBy+l );
	  HT1z_cof(k,l)=( Akx_vc(k)*HT1x_cof(k,l)+Aky_vc(l)*HT1y_cof(k,l) )/Akz_vc(k,l); 
   end;
end;

% region I

z_inc=a/NBz;
x_inc=z_inc;
zz=[0:z_inc:10*aTz-z_inc];
xx=[-aTx/2+x_inc/2:x_inc:aTx/2-x_inc/2];
zz1=zz-(10*aTz-z_inc);              % freespace1
zz2=zz;                 % freespace2

Ey_xz1=zeros(length(xx),length(zz1));
Ex_xz1=zeros(length(xx),length(zz1));
Ez_xz1=zeros(length(xx),length(zz1));
Hy_xz1=zeros(length(xx),length(zz1));
Hx_xz1=zeros(length(xx),length(zz1));
Hz_xz1=zeros(length(xx),length(zz1));

[z1 x1]=meshgrid(zz1,xx); % region 1
y=0;
  for k=1:NBx
           for l=1:NBy
           Ex_xz1=Ex_xz1+ET1x_cof(k,l)*exp(j*( Akx_vc(k)*x1+Aky_vc(l)*y-Akz_vc(k,l)*z1 ));
           Ey_xz1=Ey_xz1+ET1y_cof(k,l)*exp(j*( Akx_vc(k)*x1+Aky_vc(l)*y-Akz_vc(k,l)*z1 ));   
           Ez_xz1=Ez_xz1+ET1z_cof(k,l)*exp(j*( Akx_vc(k)*x1+Aky_vc(l)*y-Akz_vc(k,l)*z1 ));
           Hx_xz1=Hx_xz1+HT1x_cof(k,l)*exp(j*( Akx_vc(k)*x1+Aky_vc(l)*y-Akz_vc(k,l)*z1 ));
           Hy_xz1=Hy_xz1+HT1y_cof(k,l)*exp(j*( Akx_vc(k)*x1+Aky_vc(l)*y-Akz_vc(k,l)*z1 ));   
           Hz_xz1=Hz_xz1+HT1z_cof(k,l)*exp(j*( Akx_vc(k)*x1+Aky_vc(l)*y-Akz_vc(k,l)*z1 ));
           end;
        end;
       Ex_xz1=Ex_xz1+Ux*exp(j*( kx0*x1+ky0*y+kz0*z1));
       Ey_xz1=Ey_xz1+Uy*exp(j*( kx0*x1+ky0*y+kz0*z1));
       Ez_xz1=Ez_xz1+Uz*exp(j*( kx0*x1+ky0*y+kz0*z1));
       Hx_xz1=Hx_xz1+Vx*exp(j*( kx0*x1+ky0*y+kz0*z1));
       Hy_xz1=Hy_xz1+Vy*exp(j*( kx0*x1+ky0*y+kz0*z1));
       Hz_xz1=Hz_xz1+Vz*exp(j*( kx0*x1+ky0*y+kz0*z1));

% Grating region
[xlen zlen]=size(aPG_Ey_xz(:,:,1));  
aGEx_xz=zeros(length(xx),10*zlen);
aGEy_xz=zeros(length(xx),10*zlen);
aGEz_xz=zeros(length(xx),10*zlen);
aGHx_xz=zeros(length(xx),10*zlen);
aGHy_xz=zeros(length(xx),10*zlen);
aGHz_xz=zeros(length(xx),10*zlen);


   
for pbm_ind=1:pcnt
%   for pbm_ind=33:33
               
          aGEx_xz=aGEx_xz+C_p(pbm_ind)*[aPG_Ex_xz(:,:,pbm_ind)...
                aPG_Ex_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*Tz) ...
                aPG_Ex_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*2*Tz) ...
                aPG_Ex_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*3*Tz) ...
                aPG_Ex_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*4*Tz) ...
                aPG_Ex_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*5*Tz)...
                aPG_Ex_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*6*Tz) ...
                aPG_Ex_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*7*Tz) ...
                aPG_Ex_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*8*Tz) ...
                aPG_Ex_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*9*Tz) ];      
            
            aGEy_xz=aGEy_xz+C_p(pbm_ind)*[aPG_Ey_xz(:,:,pbm_ind)...
                aPG_Ey_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*Tz)...
                aPG_Ey_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*2*Tz)...
                aPG_Ey_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*3*Tz)...
                aPG_Ey_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*4*Tz)...
                aPG_Ey_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*5*Tz)...
                aPG_Ey_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*6*Tz)...
                aPG_Ey_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*7*Tz)...
                aPG_Ey_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*8*Tz)...
                aPG_Ey_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*9*Tz)]; 
      
         	aGEz_xz=aGEz_xz+C_p(pbm_ind)*[aPG_Ez_xz(:,:,pbm_ind) ...
                aPG_Ez_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*Tz)...
                aPG_Ez_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*2*Tz)...
                aPG_Ez_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*3*Tz)...
                aPG_Ez_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*4*Tz)...
                aPG_Ez_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*5*Tz)...
                aPG_Ez_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*6*Tz)...
                aPG_Ez_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*7*Tz)...
                aPG_Ez_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*8*Tz)...
                aPG_Ez_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*9*Tz)]; 
    
            aGHx_xz=aGHx_xz+C_p(pbm_ind)*[aPG_Hx_xz(:,:,pbm_ind)...
                aPG_Hx_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*Tz) ...
                aPG_Hx_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*2*Tz) ...
                aPG_Hx_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*3*Tz) ...
                aPG_Hx_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*4*Tz) ...
                aPG_Hx_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*5*Tz)...
                aPG_Hx_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*6*Tz) ...
                aPG_Hx_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*7*Tz) ...
                aPG_Hx_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*8*Tz) ...
                aPG_Hx_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*9*Tz) ];
  
            aGHy_xz=aGHy_xz+C_p(pbm_ind)*[aPG_Hy_xz(:,:,pbm_ind)...
                aPG_Hy_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*Tz)...
                aPG_Hy_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*2*Tz)...
                aPG_Hy_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*3*Tz)...
                aPG_Hy_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*4*Tz)...
                aPG_Hy_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*5*Tz)...
                aPG_Hy_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*6*Tz)...
                aPG_Hy_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*7*Tz)...
                aPG_Hy_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*8*Tz)...
                aPG_Hy_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*9*Tz)]; 
    
         	aGHz_xz=aGHz_xz+C_p(pbm_ind)*[aPG_Hz_xz(:,:,pbm_ind) ...
                aPG_Hz_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*Tz)...
                aPG_Hz_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*2*Tz)...
                aPG_Hz_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*3*Tz)...
                aPG_Hz_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*4*Tz)...
                aPG_Hz_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*5*Tz)...
                aPG_Hz_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*6*Tz)...
                aPG_Hz_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*7*Tz)...
                aPG_Hz_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*8*Tz)...
                aPG_Hz_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*9*Tz)]; 
                  
                        

    end;
            
      
figure(1); imagesc(real([Ex_xz1 aGEx_xz])); colorbar;  title('Ex-xz');
figure(2); imagesc(real([Ey_xz1 aGEy_xz])); colorbar;  title('Ey-xz');
figure(3); imagesc(real([Ez_xz1 aGEz_xz])); colorbar;  title('Ez-xz'); 

figure(4); imagesc(real([Hx_xz1 aGHx_xz])); colorbar;  title('Hx-xz'); 
figure(5); imagesc(real([Hy_xz1 aGHy_xz])); colorbar;  title('Hy-xz');
figure(6); imagesc(real([Hz_xz1 aGHz_xz])); colorbar;  title('Hz-xz');



  

