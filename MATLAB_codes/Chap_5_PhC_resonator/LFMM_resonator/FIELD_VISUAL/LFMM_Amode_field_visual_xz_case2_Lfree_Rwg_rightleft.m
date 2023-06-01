% LFMM_Amode_field_visual_xz_case2_Lfree_Rwg_rightleft

D=zeros(2*L,1);
inc_mode=33;
D(inc_mode)=1;

C_p=ARR22*D;
ET1_temp=ART21*D;
HT1_temp=-AVh*ET1_temp;

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
      
      ET1y_cof(k,l)=ET1_y((k-1)*NBy+l );
      ET1x_cof(k,l)=ET1_x((k-1)*NBy+l );
	  ET1z_cof(k,l)=( Akx_vc(k)*ET1x_cof(k,l)+Aky_vc(l)*ET1y_cof(k,l) )/Akz_vc(k,l); 
      HT1y_cof(k,l)=HT1_y((k-1)*NBy+l );
      HT1x_cof(k,l)=HT1_x((k-1)*NBy+l );
	  HT1z_cof(k,l)=( Akx_vc(k)*HT1x_cof(k,l)+Aky_vc(l)*HT1y_cof(k,l) )/Akz_vc(k,l); 
   end;
end;

SNlay=60;
% x-z Field visualization
y=0;
z_inc=a/NBz;
x_inc=z_inc;
zz=[0:z_inc:SNlay*aTz-z_inc];
xx=[-aTx/2+x_inc/2:x_inc:aTx/2-x_inc/2];
zz1=zz-(SNlay*aTz-z_inc);                      % left  freespace
zz2=zz;                                     % right freespace

% region I
Ey_xz1=zeros(length(xx),length(zz1));
Ex_xz1=zeros(length(xx),length(zz1));
Ez_xz1=zeros(length(xx),length(zz1));
Hy_xz1=zeros(length(xx),length(zz1));
Hx_xz1=zeros(length(xx),length(zz1));
Hz_xz1=zeros(length(xx),length(zz1));

% x-z Field visualization

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
        
    

% Grating region

[xlen zlen]=size(aPG_Ey_xz(:,:,1));
aGEx_xz=zeros(xlen,zlen,SNlay);
aGEy_xz=zeros(xlen,zlen,SNlay);
aGEz_xz=zeros(xlen,zlen,SNlay);
aGHx_xz=zeros(xlen,zlen,SNlay);
aGHy_xz=zeros(xlen,zlen,SNlay);
aGHz_xz=zeros(xlen,zlen,SNlay);


for laynt=1:SNlay
    for pbm_ind=1:pcnt
                        
            aGEy_xz(:,:,laynt)=aGEy_xz(:,:,laynt)+C_p(pbm_ind)*aPG_Ey_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*aTz*(laynt-1)); %Sx
            aGEx_xz(:,:,laynt)=aGEx_xz(:,:,laynt)+C_p(pbm_ind)*aPG_Ex_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*aTz*(laynt-1)); %Sx
         	aGEz_xz(:,:,laynt)=aGEz_xz(:,:,laynt)+C_p(pbm_ind)*aPG_Ez_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*aTz*(laynt-1)); 
            
            aGHy_xz(:,:,laynt)=aGHy_xz(:,:,laynt)+C_p(pbm_ind)*aPG_Hy_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*aTz*(laynt-1)); %Sx
            aGHx_xz(:,:,laynt)=aGHx_xz(:,:,laynt)+C_p(pbm_ind)*aPG_Hx_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*aTz*(laynt-1)); %Sx
         	aGHz_xz(:,:,laynt)=aGHz_xz(:,:,laynt)+C_p(pbm_ind)*aPG_Hz_xz(:,:,pbm_ind)*exp(ap_evalue(pbm_ind)*aTz*(laynt-1)); 
                  
    end;
            % incidence mode
            aGEy_xz(:,:,laynt)=aGEy_xz(:,:,laynt)+aNG_Ey_xz(:,:,inc_mode)*exp(am_evalue(inc_mode)*aTz*(laynt)); %Sx
            aGEx_xz(:,:,laynt)=aGEx_xz(:,:,laynt)+aNG_Ex_xz(:,:,inc_mode)*exp(am_evalue(inc_mode)*aTz*(laynt)); %Sx
         	aGEz_xz(:,:,laynt)=aGEz_xz(:,:,laynt)+aNG_Ez_xz(:,:,inc_mode)*exp(am_evalue(inc_mode)*aTz*(laynt)); 
            
            aGHy_xz(:,:,laynt)=aGHy_xz(:,:,laynt)+aNG_Hy_xz(:,:,inc_mode)*exp(am_evalue(inc_mode)*aTz*(laynt)); %Sx
            aGHx_xz(:,:,laynt)=aGHx_xz(:,:,laynt)+aNG_Hx_xz(:,:,inc_mode)*exp(am_evalue(inc_mode)*aTz*(laynt)); %Sx
         	aGHz_xz(:,:,laynt)=aGHz_xz(:,:,laynt)+aNG_Hz_xz(:,:,inc_mode)*exp(am_evalue(inc_mode)*aTz*(laynt)); 
           
    
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

figure(1); imagesc(real([Ex_xz1 AcaGEx_xz])); colorbar;  title('Ex-xz');
figure(2); imagesc(real([Ey_xz1 AcaGEy_xz])); colorbar;  title('Ey-xz');
figure(3); imagesc(real([Ez_xz1 AcaGEz_xz])); colorbar;  title('Ez-xz'); 

figure(4); imagesc(real([Hx_xz1 AcaGHx_xz])); colorbar;  title('Hx-xz'); 
figure(5); imagesc(real([Hy_xz1 AcaGHy_xz])); colorbar;  title('Hy-xz');
figure(6); imagesc(real([Hz_xz1 AcaGHz_xz])); colorbar;  title('Hz-xz');


