% Bmode field visualization xz plane
% freespace 4
%     |
%  grating region
%     |
% freespace 3

SNlay=60;

Ta=zeros(2*L,2*L,SNlay); % left to rignt
Ra=zeros(2*L,2*L,SNlay); % left to right

Tb=zeros(2*L,2*L,SNlay); % right to left
Rb=zeros(2*L,2*L,SNlay); % right to left

Ca=zeros(4*L,2*L,SNlay); % left to right
Cb=zeros(4*L,2*L,SNlay); % left to right


tCa=zeros(4*L,2*L,SNlay); % left to right
tCb=zeros(4*L,2*L,SNlay); % left to right


for laynt=1:SNlay
    
    Ra(:,:,laynt)=BR33;
    Ta(:,:,laynt)=BT34;
    Rb(:,:,laynt)=BR44;
    Tb(:,:,laynt)=BT43;
    
    Ca(:,:,laynt)=BCa;
    Cb(:,:,laynt)=BCb;
    
end;



   I=eye(2*L,2*L);
   T_temp1a=Ta(:,:,1);
   R_temp1a=Ra(:,:,1);
   T_temp1b=Tb(:,:,1);
   R_temp1b=Rb(:,:,1);
   
   %%% Important
   tCa=Ca;
   tCb=Cb;
   %%%

for laynt=2:SNlay
   
   %laynt
   
   T_temp2a=Ta(:,:,laynt);
   R_temp2a=Ra(:,:,laynt);
   T_temp2b=Tb(:,:,laynt);
   R_temp2b=Rb(:,:,laynt);

	RRa=(R_temp1a+T_temp1b*inv(I-R_temp2a*R_temp1b)*R_temp2a*T_temp1a);
	TTa=T_temp2a*inv(I-R_temp1b*R_temp2a)*T_temp1a;

	RRb=(R_temp2b+T_temp2a*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b);
   TTb=T_temp1b*inv(I-R_temp2a*R_temp1b)*T_temp2b;
   
   
   for k=1:laynt-1
   	
   tCa(:,:,k)=Ca(:,:,k)+Cb(:,:,k)*inv(I-R_temp2a*R_temp1b)*R_temp2a*T_temp1a;
   tCb(:,:,k)=Cb(:,:,k)*inv(I-R_temp2a*R_temp1b)*T_temp2b;

   end; % for k
   
   tCa(:,:,laynt)=Ca(:,:,laynt)*inv(I-R_temp1b*R_temp2a)*T_temp1a;
   tCb(:,:,laynt)=Cb(:,:,laynt)+Ca(:,:,laynt)*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b;
   
   T_temp1a=TTa;
   R_temp1a=RRa;
   T_temp1b=TTb;
   R_temp1b=RRb;
   
   Ca=tCa;
   Cb=tCb; 
end; % laynt

   SBT34=T_temp1a;
   SBR33=R_temp1a;
   SBT43=T_temp1b;
   SBR44=R_temp1b;
   SBCa=Ca;
   SBCb=Cb;

% Field visualization 

x_inc=a/9;
z_inc=x_inc;
xx=[0:x_inc:10*bTx-x_inc];
zz=[-bTz/2+z_inc/2:z_inc:bTz/2-z_inc/2];
xx3=xx-(10*bTx-x_inc);              % freespace1
xx4=xx;                 % freespace2

% region 3
Ey_xz3=zeros(length(xx3),length(zz));
Ex_xz3=zeros(length(xx3),length(zz));
Ez_xz3=zeros(length(xx3),length(zz));
% region 4
Ey_xz4=zeros(length(xx4),length(zz));
Ex_xz4=zeros(length(xx4),length(zz));
Ez_xz4=zeros(length(xx4),length(zz));

y=0;
Uy=1;
Ux=0;
D=zeros(2*L,1);
d11=zeros(L,1);
d12=zeros(L,1);
centx=nx*NBy+ny+1;
centy=nx*NBy+ny+1;
d11(centx)=Uy;  % Kx of incident beam
d12(centy)=Ux;
D(1:L)=d11;
D(L+1:2*L)=d12;

T3_temp=SBR33*D;
T4_temp=SBT34*D;

for laynt=1:SNlay

C_temp(:,laynt)=Ca(:,:,laynt)*D;
C_p(:,laynt)=C_temp(1:2*L,laynt);          % Positive coupling coefficients
C_n(:,laynt)=C_temp(2*L+1:4*L,laynt);      % Negative coupling coefficients

end;

T3_y=zeros(L,1);
T3_z=zeros(L,1);
T3_y=T3_temp(1:L);
T3_z=T3_temp(L+1:2*L);

T4_z=zeros(L,1);
T4_y=zeros(L,1);
T4_y=T4_temp(1:L);
T4_z=T4_temp(L+1:2*L);
% 
% T3x_cof=zeros(NBz,NBy);
% T3y_cof=zeros(NBz,NBy);
% T3z_cof=zeros(NBz,NBy);
% 
% T4x_cof=zeros(NBz,NBy);
% T4y_cof=zeros(NBz,NBy);
% T4z_cof=zeros(NBz,NBy);

T3x_cof=zeros(NBx,NBy);
T3y_cof=zeros(NBx,NBy);
T3z_cof=zeros(NBx,NBy);

T4x_cof=zeros(NBx,NBy);
T4y_cof=zeros(NBx,NBy);
T4z_cof=zeros(NBx,NBy);


for k=1:NBx
   for l=1:NBy
      
      T4y_cof(k,l)=T4_y((k-1)*NBy+l);
      T4z_cof(k,l)=T4_z((k-1)*NBy+l );
      T4x_cof(k,l)=-( Bkz_vc(k)*T4z_cof(k,l)+Bky_vc(l)*T4y_cof(k,l) )/Bkx_vc(k,l); 

      T3y_cof(k,l)=T3_y((k-1)*NBy+l );
      T3z_cof(k,l)=T3_z((k-1)*NBy+l );
	  T3x_cof(k,l)=( Bkz_vc(k)*T3z_cof(k,l)+Bky_vc(l)*T3y_cof(k,l) )/Bkx_vc(k,l); 
   end;
end;


% diffraction efficiency

DEt4=zeros(NBx,NBy);
DEt3=zeros(NBx,NBy);
for k=1:NBx
   for l=1:NBy
      
      DEt3(k,l)=abs(abs(T3x_cof(k,l))^2+abs(T3y_cof(k,l))^2+abs(T3z_cof(k,l))^2)*real(Bkx_vc(k,l)/(kz0));
      DEt4(k,l)=abs(abs(T4x_cof(k,l))^2+abs(T4y_cof(k,l))^2+abs(T4z_cof(k,l))^2)*real(Bkx_vc(k,l)/(kz0));
      
   end;   
end;
   

total_energy=sum(sum(DEt3))+sum(sum(DEt4));


% x-z Field visualization

% region 3
[z3 x3]=meshgrid(zz,xx3); % region 3

  for k=1:NBx
           for l=1:NBy
           Ey_xz3=Ey_xz3+T3y_cof(k,l)*exp(j*( Bkz_vc(k)*z3+Bky_vc(l)*y-Bkx_vc(k,l)*x3));   
           Ex_xz3=Ex_xz3+T3x_cof(k,l)*exp(j*( Bkz_vc(k)*z3+Bky_vc(l)*y-Bkx_vc(k,l)*x3 ));
           Ez_xz3=Ez_xz3+T3z_cof(k,l)*exp(j*( Bkz_vc(k)*z3+Bky_vc(l)*y-Bkx_vc(k,l)*x3 ));
           end;
        end;
       Ey_xz3=Ey_xz3+Uy*exp(j*( kx0*z3+ky0*y+kz0*x3));
       Ex_xz3=Ex_xz3+Uz*exp(j*( kx0*z3+ky0*y+kz0*x3));
       Ez_xz3=Ez_xz3+Ux*exp(j*( kx0*z3+ky0*y+kz0*x3));
  
       
% region 4
[z4 x4]=meshgrid(zz,xx4); % region 4
    for k=1:NBx
         for l=1:NBy
           Ey_xz4=Ey_xz4+T4y_cof(k,l)*exp(j*( Bkz_vc(k)*z4+Bky_vc(l)*y+Bkx_vc(k,l)*x4 ));   
           Ex_xz4=Ex_xz4+T4x_cof(k,l)*exp(j*( Bkz_vc(k)*z4+Bky_vc(l)*y+Bkx_vc(k,l)*x4 ));
           Ez_xz4=Ez_xz4+T4z_cof(k,l)*exp(j*( Bkz_vc(k)*z4+Bky_vc(l)*y+Bkx_vc(k,l)*x4 ));
           end;
        end;
 
% Grating region
xx=[0:x_inc:bTx-x_inc];
zz=[-bTz/2+z_inc/2:z_inc:bTz/2-z_inc/2];

bGEx_xz=zeros(length(xx),length(zz),SNlay);
bGEy_xz=zeros(length(xx),length(zz),SNlay);
bGEz_xz=zeros(length(xx),length(zz),SNlay);

for laynt=1:SNlay
   for pbm_ind=1:pcnt
                        
            bGEy_xz(:,:,laynt)=bGEy_xz(:,:,laynt)+C_p(pbm_ind,laynt)*bPG_Ey_xz(:,:,pbm_ind); %Sx
            bGEx_xz(:,:,laynt)=bGEx_xz(:,:,laynt)+C_p(pbm_ind,laynt)*bPG_Ex_xz(:,:,pbm_ind); %Sx
         	bGEz_xz(:,:,laynt)=bGEz_xz(:,:,laynt)+C_p(pbm_ind,laynt)*bPG_Ez_xz(:,:,pbm_ind); 
                  
    end;
            
      
   for pbm_ind=1:mcnt
                        
            bGEy_xz(:,:,laynt)=bGEy_xz(:,:,laynt)+C_n(pbm_ind,laynt)*bNG_Ey_xz(:,:,pbm_ind); %Sx
            bGEx_xz(:,:,laynt)=bGEx_xz(:,:,laynt)+C_n(pbm_ind,laynt)*bNG_Ex_xz(:,:,pbm_ind); %Sx
         	bGEz_xz(:,:,laynt)=bGEz_xz(:,:,laynt)+C_n(pbm_ind,laynt)*bNG_Ez_xz(:,:,pbm_ind); 
                  
    end;
end;
    
BcaGEx_xz=bGEx_xz(:,:,1);
BcaGEy_xz=bGEy_xz(:,:,1);
BcaGEz_xz=bGEz_xz(:,:,1);

for laynt=2:SNlay
    BcaGEx_xz=[BcaGEx_xz; bGEx_xz(:,:,laynt)];
    BcaGEy_xz=[BcaGEy_xz; bGEy_xz(:,:,laynt)];
    BcaGEz_xz=[BcaGEz_xz; bGEz_xz(:,:,laynt)];
end;

figure(1); imagesc(real([Ex_xz3; BcaGEx_xz; Ex_xz4])); set(gca,'ydir','normal'); title('Ex-xz');  colorbar;
figure(2); imagesc(real([Ey_xz3; BcaGEy_xz; Ey_xz4])); set(gca,'ydir','normal'); title('Ey-xz'); colorbar;
figure(3); imagesc(real([Ez_xz3; BcaGEz_xz; Ez_xz4])); set(gca,'ydir','normal'); title('Ez-xz'); colorbar;

% for kk=0:0.2:100*2*pi
% figure(4); imagesc((real([Ey_xz3; BcaGEy_xz; Ey_xz4]*(exp(-i*kk)))));  caxis([-2 2]);
% pause(0.1);
% end;




