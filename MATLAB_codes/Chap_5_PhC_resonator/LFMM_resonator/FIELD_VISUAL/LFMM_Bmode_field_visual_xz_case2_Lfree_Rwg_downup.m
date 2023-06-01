% LFMM_Bmode_field_visual_xz_case2_Lfree_Rwg_downup

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
Hinc=BZh*Einc;
Vy=Hinc(centx);
Vz=Hinc(L+centx);
Vx=-(kx0*Vx+ky0*Vy)/kz0;

ET3_temp=BRR33*Einc;
HT3_temp=-BZh*ET3_temp;
C_p=BRT34*Einc;

ET3_y=zeros(L,1);
ET3_z=zeros(L,1);
HT3_y=zeros(L,1);
HT3_z=zeros(L,1);

ET3_y=ET3_temp(1:L);
ET3_z=ET3_temp(L+1:2*L);
HT3_y=HT3_temp(1:L);
HT3_z=HT3_temp(L+1:2*L);

ET3x_cof=zeros(NBx,NBy);
ET3y_cof=zeros(NBx,NBy);
ET3z_cof=zeros(NBx,NBy);
HT3x_cof=zeros(NBx,NBy);
HT3y_cof=zeros(NBx,NBy);
HT3z_cof=zeros(NBx,NBy);

for k=1:NBx
   for l=1:NBy
      
      ET3y_cof(k,l)=ET3_y((k-1)*NBy+l );
      ET3z_cof(k,l)=ET3_z((k-1)*NBy+l );
	  ET3x_cof(k,l)=( Bkz_vc(k)*ET3z_cof(k,l)+Bky_vc(l)*ET3y_cof(k,l) )/Bkx_vc(k,l); 
      HT3y_cof(k,l)=HT3_y((k-1)*NBy+l );
      HT3z_cof(k,l)=HT3_z((k-1)*NBy+l );
	  HT3x_cof(k,l)=( Bkz_vc(k)*HT3z_cof(k,l)+Bky_vc(l)*HT3y_cof(k,l) )/Bkx_vc(k,l); 
   
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


% x-z Field visualization
y=0;

x_inc=a/NBx;
z_inc=x_inc;
xx=[0:x_inc:SNlay*bTx-x_inc];
zz=[-bTz/2+z_inc/2:z_inc:bTz/2-z_inc/2];
xx3=xx-(SNlay*bTx-x_inc);                      % down  freespace
xx4=xx;                                     % up freespace

% region 3
Ey_xz3=zeros(length(xx3),length(zz));
Ex_xz3=zeros(length(xx3),length(zz));
Ez_xz3=zeros(length(xx3),length(zz));
Hy_xz3=zeros(length(xx3),length(zz));
Hx_xz3=zeros(length(xx3),length(zz));
Hz_xz3=zeros(length(xx3),length(zz));


% region 3
[z3 x3]=meshgrid(zz,xx3); 
  for k=1:NBx
           for l=1:NBy
           Ey_xz3=Ey_xz3+ET3y_cof(k,l)*exp(j*( Bkz_vc(k)*z3+Bky_vc(l)*y-Bkx_vc(k,l)*x3 ));   
           Ex_xz3=Ex_xz3+ET3x_cof(k,l)*exp(j*( Bkz_vc(k)*z3+Bky_vc(l)*y-Bkx_vc(k,l)*x3 ));
           Ez_xz3=Ez_xz3+ET3z_cof(k,l)*exp(j*( Bkz_vc(k)*z3+Bky_vc(l)*y-Bkx_vc(k,l)*x3 ));
           
           Hy_xz3=Hy_xz3+HT3y_cof(k,l)*exp(j*( Bkz_vc(k)*z3+Bky_vc(l)*y-Bkx_vc(k,l)*x3 ));   
           Hx_xz3=Hx_xz3+HT3x_cof(k,l)*exp(j*( Bkz_vc(k)*z3+Bky_vc(l)*y-Bkx_vc(k,l)*x3 ));
           Hz_xz3=Hz_xz3+HT3z_cof(k,l)*exp(j*( Bkz_vc(k)*z3+Bky_vc(l)*y-Bkx_vc(k,l)*x3 ));
           end;
        end;
        
       Ey_xz3=Ey_xz3+Uy*exp(j*( kx0*z3+ky0*y+kz0*x3));
       Ex_xz3=Ex_xz3+Ux*exp(j*( kx0*z3+ky0*y+kz0*x3));
       Ez_xz3=Ez_xz3+Uz*exp(j*( kx0*z3+ky0*y+kz0*x3));
       Hy_xz3=Hy_xz3+Vy*exp(j*( kx0*z3+ky0*y+kz0*x3));
       Hx_xz3=Hx_xz3+Vx*exp(j*( kx0*z3+ky0*y+kz0*x3));
       Hz_xz3=Hz_xz3+Vz*exp(j*( kx0*z3+ky0*y+kz0*x3));



% Grating region

SNlay=60;

[xlen zlen]=size(bPG_Ey_xz(:,:,1));
bGEx_xz=zeros(xlen,zlen,SNlay);
bGEy_xz=zeros(xlen,zlen,SNlay);
bGEz_xz=zeros(xlen,zlen,SNlay);
bGHx_xz=zeros(xlen,zlen,SNlay);
bGHy_xz=zeros(xlen,zlen,SNlay);
bGHz_xz=zeros(xlen,zlen,SNlay);

for laynt=1:SNlay
   for pbm_ind=1:pcnt
                        
            bGEy_xz(:,:,laynt)=bGEy_xz(:,:,laynt)+C_p(pbm_ind)*bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx*(laynt-1)); %Sx
            bGEx_xz(:,:,laynt)=bGEx_xz(:,:,laynt)+C_p(pbm_ind)*bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx*(laynt-1)); %Sx
         	bGEz_xz(:,:,laynt)=bGEz_xz(:,:,laynt)+C_p(pbm_ind)*bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx*(laynt-1)); 
            
            bGHy_xz(:,:,laynt)=bGHy_xz(:,:,laynt)+C_p(pbm_ind)*bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx*(laynt-1)); %Sx
            bGHx_xz(:,:,laynt)=bGHx_xz(:,:,laynt)+C_p(pbm_ind)*bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx*(laynt-1)); %Sx
         	bGHz_xz(:,:,laynt)=bGHz_xz(:,:,laynt)+C_p(pbm_ind)*bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx*(laynt-1)); 
                  
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

figure(1); imagesc(real([Ex_xz3; BcaGEx_xz;])); title('Ex-xz'); colorbar;
figure(2); imagesc(real([Ey_xz3; BcaGEy_xz;])); title('Ey-xz'); colorbar;
figure(3); imagesc(real([Ez_xz3; BcaGEz_xz;])); title('Ez-xz'); colorbar;
figure(4); imagesc(real([Hx_xz3; BcaGHx_xz;])); title('Hx-xz'); colorbar;
figure(5); imagesc(real([Hy_xz3; BcaGHy_xz;])); title('Hy-xz'); colorbar;
figure(6); imagesc(real([Hz_xz3; BcaGHz_xz;])); title('Hz-xz'); colorbar;

% 
% for kk=0:0.2:10*2*pi
% figure(7); imagesc((real([Ey_xz1 AcaGEy_xz Ey_xz2]*(exp(-i*kk)))));  caxis([-3 3]);
% end;



