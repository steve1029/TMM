% LFMM_Bmode_field_visual_xz_case3_Lwg_Rfree_updown

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
Hinc=-BZh*Einc;
Vy=Hinc(centx);
Vz=Hinc(L+centx);
Vx=(kx0*Vx+ky0*Vy)/kz0;

ET4_temp=BLR44*Einc;
HT4_temp=BZh*ET4_temp;
C_n=BLT43*Einc;

ET4_y=zeros(L,1);
ET4_z=zeros(L,1);
HT4_y=zeros(L,1);
HT4_z=zeros(L,1);

ET4_y=ET4_temp(1:L);
ET4_z=ET4_temp(L+1:2*L);
HT4_y=HT4_temp(1:L);
HT4_z=HT4_temp(L+1:2*L);

ET4x_cof=zeros(NBx,NBy);
ET4y_cof=zeros(NBx,NBy);
ET4z_cof=zeros(NBx,NBy);
HT4x_cof=zeros(NBx,NBy);
HT4y_cof=zeros(NBx,NBy);
HT4z_cof=zeros(NBx,NBy);

for k=1:NBx
   for l=1:NBy
      
      ET4y_cof(k,l)=ET4_y((k-1)*NBy+l );
      ET4z_cof(k,l)=ET4_z((k-1)*NBy+l );
	  ET4x_cof(k,l)=-( Bkz_vc(k)*ET4z_cof(k,l)+Bky_vc(l)*ET4y_cof(k,l) )/Bkx_vc(k,l); 
      HT4y_cof(k,l)=HT4_y((k-1)*NBy+l );
      HT4z_cof(k,l)=HT4_z((k-1)*NBy+l );
	  HT4x_cof(k,l)=-( Bkz_vc(k)*HT4z_cof(k,l)+Bky_vc(l)*HT4y_cof(k,l) )/Bkx_vc(k,l); 
   
   end;
end;

% % diffraction efficiency
% 
% DEt4=zeros(NBx,NBy);
% DEt3=zeros(NBx,NBy);
% for k=1:NBx
%    for l=1:NBy
%       
%       DEt3(k,l)=abs(abs(ET3x_cof(k,l))^2+abs(ET3y_cof(k,l))^2+abs(ET3z_cof(k,l))^2)*real(Bkx_vc(k,l)/(kz0));
%       DEt4(k,l)=abs(abs(ET4x_cof(k,l))^2+abs(ET4y_cof(k,l))^2+abs(ET4z_cof(k,l))^2)*real(Bkx_vc(k,l)/(kz0));
%       
%    end;   
% end;
%   
% total_energy=sum(sum(DEt3))+sum(sum(DEt4));
% 

SNlay=60;
% x-z Field visualization
y=0;

x_inc=a/NBx;
z_inc=x_inc;
xx=[0:x_inc:SNlay*bTx-x_inc];
zz=[-bTz/2+z_inc/2:z_inc:bTz/2-z_inc/2];
xx3=xx-(SNlay*bTx-x_inc);                       % down  freespace
xx4=xx;                                         % up freespace

% region 4
Ey_xz4=zeros(length(xx4),length(zz));
Ex_xz4=zeros(length(xx4),length(zz));
Ez_xz4=zeros(length(xx4),length(zz));
Hy_xz4=zeros(length(xx4),length(zz));
Hx_xz4=zeros(length(xx4),length(zz));
Hz_xz4=zeros(length(xx4),length(zz));


% region 4
[z4 x4]=meshgrid(zz,xx4); 
  for k=1:NBx
           for l=1:NBy
           Ey_xz4=Ey_xz4+ET4y_cof(k,l)*exp(j*( Bkz_vc(k)*z4+Bky_vc(l)*y+Bkx_vc(k,l)*x4 ));   
           Ex_xz4=Ex_xz4+ET4x_cof(k,l)*exp(j*( Bkz_vc(k)*z4+Bky_vc(l)*y+Bkx_vc(k,l)*x4 ));
           Ez_xz4=Ez_xz4+ET4z_cof(k,l)*exp(j*( Bkz_vc(k)*z4+Bky_vc(l)*y+Bkx_vc(k,l)*x4 ));
           
           Hy_xz4=Hy_xz4+HT4y_cof(k,l)*exp(j*( Bkz_vc(k)*z4+Bky_vc(l)*y+Bkx_vc(k,l)*x4 ));   
           Hx_xz4=Hx_xz4+HT4x_cof(k,l)*exp(j*( Bkz_vc(k)*z4+Bky_vc(l)*y+Bkx_vc(k,l)*x4 ));
           Hz_xz4=Hz_xz4+HT4z_cof(k,l)*exp(j*( Bkz_vc(k)*z4+Bky_vc(l)*y+Bkx_vc(k,l)*x4 ));
           end;
        end;
        
       Ey_xz4=Ey_xz4+Uy*exp(j*( kx0*z4+ky0*y-kz0*x4));
       Ex_xz4=Ex_xz4+Ux*exp(j*( kx0*z4+ky0*y-kz0*x4));
       Ez_xz4=Ez_xz4+Uz*exp(j*( kx0*z4+ky0*y-kz0*x4));
       Hy_xz4=Hy_xz4+Vy*exp(j*( kx0*z4+ky0*y-kz0*x4));
       Hx_xz4=Hx_xz4+Vx*exp(j*( kx0*z4+ky0*y-kz0*x4));
       Hz_xz4=Hz_xz4+Vz*exp(j*( kx0*z4+ky0*y-kz0*x4));


% Grating region
[xlen zlen]=size(bNG_Ey_xz(:,:,1));
bGEx_xz=zeros(xlen,zlen,SNlay);
bGEy_xz=zeros(xlen,zlen,SNlay);
bGEz_xz=zeros(xlen,zlen,SNlay);
bGHx_xz=zeros(xlen,zlen,SNlay);
bGHy_xz=zeros(xlen,zlen,SNlay);
bGHz_xz=zeros(xlen,zlen,SNlay);

for laynt=1:SNlay
   for pbm_ind=1:mcnt
                        
            bGEy_xz(:,:,laynt)=bGEy_xz(:,:,laynt)+C_n(pbm_ind)*bNG_Ey_xz(:,:,pbm_ind)*exp(-bm_evalue(pbm_ind)*bTx*(SNlay-laynt)); %Sx
            bGEx_xz(:,:,laynt)=bGEx_xz(:,:,laynt)+C_n(pbm_ind)*bNG_Ex_xz(:,:,pbm_ind)*exp(-bm_evalue(pbm_ind)*bTx*(SNlay-laynt)); %Sx
         	bGEz_xz(:,:,laynt)=bGEz_xz(:,:,laynt)+C_n(pbm_ind)*bNG_Ez_xz(:,:,pbm_ind)*exp(-bm_evalue(pbm_ind)*bTx*(SNlay-laynt)); 
            
            bGHy_xz(:,:,laynt)=bGHy_xz(:,:,laynt)+C_n(pbm_ind)*bNG_Hy_xz(:,:,pbm_ind)*exp(-bm_evalue(pbm_ind)*bTx*(SNlay-laynt)); %Sx
            bGHx_xz(:,:,laynt)=bGHx_xz(:,:,laynt)+C_n(pbm_ind)*bNG_Hx_xz(:,:,pbm_ind)*exp(-bm_evalue(pbm_ind)*bTx*(SNlay-laynt)); %Sx
         	bGHz_xz(:,:,laynt)=bGHz_xz(:,:,laynt)+C_n(pbm_ind)*bNG_Hz_xz(:,:,pbm_ind)*exp(-bm_evalue(pbm_ind)*bTx*(SNlay-laynt)); 
                  
    end;
            
      

end; % for laynt


BcaGEx_xz=bGEx_xz(:,:,1);
BcaGEy_xz=bGEy_xz(:,:,1);
BcaGEz_xz=bGEz_xz(:,:,1);

BcaGHx_xz=bGHx_xz(:,:,1);
BcaGHy_xz=bGHy_xz(:,:,1);
BcaGHz_xz=bGHz_xz(:,:,1);

for laynt=2:SNlay
    BcaGEx_xz=[BcaGEx_xz; bGEx_xz(:,:,laynt)];
    BcaGEy_xz=[BcaGEy_xz; bGEy_xz(:,:,laynt)];
    BcaGEz_xz=[BcaGEz_xz; bGEz_xz(:,:,laynt)];
    
    BcaGHx_xz=[BcaGHx_xz; bGHx_xz(:,:,laynt)];
    BcaGHy_xz=[BcaGHy_xz; bGHy_xz(:,:,laynt)];
    BcaGHz_xz=[BcaGHz_xz; bGHz_xz(:,:,laynt)];
end;

figure(1); imagesc(real([BcaGEx_xz; Ex_xz4])); title('Ex-xz'); colorbar;
figure(2); imagesc(real([BcaGEy_xz; Ey_xz4])); title('Ey-xz'); colorbar;
figure(3); imagesc(real([BcaGEz_xz; Ez_xz4])); title('Ez-xz'); colorbar;
figure(4); imagesc(real([BcaGHx_xz; Hx_xz4])); title('Hx-xz'); colorbar;
figure(5); imagesc(real([BcaGHy_xz; Hy_xz4])); title('Hy-xz'); colorbar;
figure(6); imagesc(real([BcaGHz_xz; Hz_xz4])); title('Hz-xz'); colorbar;

% 
% for kk=0:0.2:10*2*pi
% figure(7); imagesc((real([Ey_xz1 AcaGEy_xz Ey_xz2]*(exp(-i*kk)))));  caxis([-3 3]);
% end;



