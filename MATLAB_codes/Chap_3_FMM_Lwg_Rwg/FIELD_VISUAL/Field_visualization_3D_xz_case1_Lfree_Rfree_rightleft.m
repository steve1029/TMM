%% Field visualization case 1 right-left

xx=[-Tx/2:Tx*0.005:Tx/2]';
z_start=-2*um;
z_end=5*um;
z_inc=20*nano;
zz=(z_start:z_inc:z_end);
zn=length(zz);

zrev=(z_end:-z_inc:z_start);


% x-z plane
Ex_xz=zeros(length(xx),length(zz));
Ey_xz=zeros(length(xx),length(zz));
Ez_xz=zeros(length(xx),length(zz));
Hx_xz=zeros(length(xx),length(zz));
Hy_xz=zeros(length(xx),length(zz));
Hz_xz=zeros(length(xx),length(zz));

rg1cnt=0;               % region 1
rgGcnt=zeros(Nlay,1);   % grating region 
rg2cnt=0;               % region 2

for znt=1:zn
    
    zv=zz(znt);
    
    if zv <= 0
      rg1cnt=rg1cnt+1;
    end;
    
    if (0 < zv) & (zv <= ac_thick(Nlay))
            
       for laynt=1:Nlay
             if  ( ( ac_thick(laynt)-lay_thick(laynt) ) < zv ) & ( zv <= ac_thick(laynt)  )
              rgGcnt(laynt)=rgGcnt(laynt)+1; 
             end;
       end;
    
    end;
    
     if ( ac_thick(Nlay) < zv )
        rg2cnt=rg2cnt+1;
    end;
    
    
end; % for znt

 % region 1
 zz1ind(1,1)=1;
 zz1ind(1,2)=rg1cnt;
 
 % grating region
 for laynt=1:Nlay
     
     if laynt==1
 zzGind(laynt,1)=rg1cnt+1;
 zzGind(laynt,2)=zzGind(laynt,1)+rgGcnt(laynt)-1 ;
 
     else
 zzGind(laynt,1)=zzGind(laynt-1,2)+1;
 zzGind(laynt,2)=zzGind(laynt,1)+rgGcnt(laynt)-1 ;
     end;
 end;
 
 fc_ind=zzGind-rg1cnt; % field calculation index for grating region
 
 % region 2
 zz2ind(1,1)=rg1cnt+sum(rgGcnt)+1;
 zz2ind(1,2)=rg1cnt+sum(rgGcnt)+rg2cnt;
 
 % region 1
  Ey_xz1=zeros(length(xx),rg1cnt);
  Ex_xz1=zeros(length(xx),rg1cnt);
  Ez_xz1=zeros(length(xx),rg1cnt);
  Hy_xz1=zeros(length(xx),rg1cnt);
  Hx_xz1=zeros(length(xx),rg1cnt);
  Hz_xz1=zeros(length(xx),rg1cnt);     

 % Grating region
 
  GEy_xz=zeros(length(xx),sum(rgGcnt));
  GEx_xz=zeros(length(xx),sum(rgGcnt));
  GEz_xz=zeros(length(xx),sum(rgGcnt));
  GHy_xz=zeros(length(xx),sum(rgGcnt));
  GHx_xz=zeros(length(xx),sum(rgGcnt));
  GHz_xz=zeros(length(xx),sum(rgGcnt));     

  % region 2
  
  Ey_xz2=zeros(length(xx),rg2cnt);
  Ex_xz2=zeros(length(xx),rg2cnt);
  Ez_xz2=zeros(length(xx),rg2cnt);
  Hy_xz2=zeros(length(xx),rg2cnt);
  Hx_xz2=zeros(length(xx),rg2cnt);
  Hz_xz2=zeros(length(xx),rg2cnt);     

  Iden=zeros(2*L,2*L);
for k=1:2*L
Iden(k,k)=1;   
end;
       
%% .    input field - plane wave


%input_pol_mode='TM';
input_pol_mode='TE';

if input_pol_mode == 'TM'

Uy=tm_Uy;
Ux=tm_Ux;
Uz=tm_Uz;

Einc=zeros(2*L,1);
d11=zeros(L,1);
d12=zeros(L,1);
centx=nx*NBy+ny+1;
centy=nx*NBy+ny+1;

d11(centy)=Uy;  % Kx of incident beam
d12(centx)=Ux;

Einc(1:L)=d11;
Einc(L+1:2*L)=d12;
Hinc=-VhII*Einc;
Vy=Hinc(centx);
Vx=Hinc(L+centx);
Vz=(kix*Vx+kiy*Vy)/kfz;


tm_Vy=Vy;
tm_Vx=Vx;
tm_Vz=Vz;

end;

if input_pol_mode == 'TE'

Uy=te_Uy;
Ux=te_Ux;
Uz=te_Uz;

Einc=zeros(2*L,1);
d11=zeros(L,1);
d12=zeros(L,1);
centx=nx*NBy+ny+1;
centy=nx*NBy+ny+1;

d11(centy)=Uy;  % Kx of incident beam
d12(centx)=Ux;

Einc(1:L)=d11;
Einc(L+1:2*L)=d12;
Hinc=-VhII*Einc;
Vy=Hinc(centx);
Vx=Hinc(L+centx);
Vz=(kix*Vx+kiy*Vy)/kfz;

te_Vy=Vy;
te_Vx=Vx;
te_Vz=Vz;

end;

ET1_temp=TTb*Einc;
HT1_temp=-VhI*ET1_temp;
ET2_temp=RRb*Einc;
HT2_temp=VhII*ET2_temp;


for laynt=1:Nlay
C_temp(:,laynt)=Cb(:,:,laynt)*Einc;
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
      ET2z_cof(k,l)=-( kx_tra(k)*ET2x_cof(k,l)+ky_tra(l)*ET2y_cof(k,l) )/kz_tra(k,l); 
      HT2y_cof(k,l)=HT2_y((k-1)*NBy+l);
      HT2x_cof(k,l)=HT2_x((k-1)*NBy+l );
      HT2z_cof(k,l)=-( kx_tra(k)*HT2x_cof(k,l)+ky_tra(l)*HT2y_cof(k,l) )/kz_tra(k,l); 

      ET1y_cof(k,l)=ET1_y((k-1)*NBy+l );
      ET1x_cof(k,l)=ET1_x((k-1)*NBy+l );
	  ET1z_cof(k,l)=( kx_ref(k)*ET1x_cof(k,l)+ky_ref(l)*ET1y_cof(k,l) )/kz_ref(k,l); 
      HT1y_cof(k,l)=HT1_y((k-1)*NBy+l );
      HT1x_cof(k,l)=HT1_x((k-1)*NBy+l );
	  HT1z_cof(k,l)=( kx_ref(k)*HT1x_cof(k,l)+ky_ref(l)*HT1y_cof(k,l) )/kz_ref(k,l); 
   
   end;
end;


% diffraction efficiency

DEt1=zeros(NBx,NBy);
DEt2=zeros(NBx,NBy);
for k=1:NBx
   for l=1:NBy
      
      DEt1(k,l)=abs(abs(ET1x_cof(k,l))^2+abs(ET1y_cof(k,l))^2+abs(ET1z_cof(k,l))^2)*real(kz_ref(k,l)/(kfz));
      DEt2(k,l)=abs(abs(ET2x_cof(k,l))^2+abs(ET2y_cof(k,l))^2+abs(ET2z_cof(k,l))^2)*real(kz_tra(k,l)/(kfz));
      
   end;   
end;
   
total_energy=sum(sum(DEt1))+sum(sum(DEt2));

% x-z Field visualization
y=0;
% region 1
[z1 x1]=meshgrid(zz(zz1ind(1):zz1ind(2)),xx); % region 1
  for k=1:NBx
           for l=1:NBy
           Ey_xz1=Ey_xz1+ET1y_cof(k,l)*exp(j*( kx_ref(k)*x1+ky_ref(l)*y-kz_ref(k,l)*z1 ));   
           Ex_xz1=Ex_xz1+ET1x_cof(k,l)*exp(j*( kx_ref(k)*x1+ky_ref(l)*y-kz_ref(k,l)*z1 ));
           Ez_xz1=Ez_xz1+ET1z_cof(k,l)*exp(j*( kx_ref(k)*x1+ky_ref(l)*y-kz_ref(k,l)*z1 ));
           
           Hy_xz1=Hy_xz1+HT1y_cof(k,l)*exp(j*( kx_ref(k)*x1+ky_ref(l)*y-kz_ref(k,l)*z1 ));   
           Hx_xz1=Hx_xz1+HT1x_cof(k,l)*exp(j*( kx_ref(k)*x1+ky_ref(l)*y-kz_ref(k,l)*z1 ));
           Hz_xz1=Hz_xz1+HT1z_cof(k,l)*exp(j*( kx_ref(k)*x1+ky_ref(l)*y-kz_ref(k,l)*z1 ));
           end;
        end;
        
%    Ey_xz1=Ey_xz1+Uy*exp(j*( kix*x1+kiy*y+kiz*z1));
%        Ex_xz1=Ex_xz1+Ux*exp(j*( kix*x1+kiy*y+kiz*z1));
%        Ez_xz1=Ez_xz1+Uz*exp(j*( kix*x1+kiy*y+kiz*z1));
%        Hy_xz1=Hy_xz1+Vy*exp(j*( kix*x1+kiy*y+kiz*z1));
%        Hx_xz1=Hx_xz1+Vx*exp(j*( kix*x1+kiy*y+kiz*z1));
%        Hz_xz1=Hz_xz1+Vz*exp(j*( kix*x1+kiy*y+kiz*z1));



% region 3

[z2 x2]=meshgrid(zz(zz2ind(1):zz2ind(2)),xx); % region 2
    for k=1:NBx
         for l=1:NBy
           Ey_xz2=Ey_xz2+ET2y_cof(k,l)*exp(j*kz_tra(k,l)*(z2-ac_thick(Nlay))).*exp(j*( kx_tra(k)*x2+ky_tra(l)*y ));   
           Ex_xz2=Ex_xz2+ET2x_cof(k,l)*exp(j*kz_tra(k,l)*(z2-ac_thick(Nlay))).*exp(j*( kx_tra(k)*x2+ky_tra(l)*y ));
           Ez_xz2=Ez_xz2+ET2z_cof(k,l)*exp(j*kz_tra(k,l)*(z2-ac_thick(Nlay))).*exp(j*( kx_tra(k)*x2+ky_tra(l)*y ));
           Hy_xz2=Hy_xz2+HT2y_cof(k,l)*exp(j*kz_tra(k,l)*(z2-ac_thick(Nlay))).*exp(j*( kx_tra(k)*x2+ky_tra(l)*y ));   
           Hx_xz2=Hx_xz2+HT2x_cof(k,l)*exp(j*kz_tra(k,l)*(z2-ac_thick(Nlay))).*exp(j*( kx_tra(k)*x2+ky_tra(l)*y ));
           Hz_xz2=Hz_xz2+HT2z_cof(k,l)*exp(j*kz_tra(k,l)*(z2-ac_thick(Nlay))).*exp(j*( kx_tra(k)*x2+ky_tra(l)*y ));

            end;
        end;
        
        Ey_xz2=Ey_xz2+Uy*exp(j*( kix*x2+kiy*y-kfz*(z2-ac_thick(Nlay))));
        Ex_xz2=Ex_xz2+Ux*exp(j*( kix*x2+kiy*y-kfz*(z2-ac_thick(Nlay))));
        Ez_xz2=Ez_xz2+Uz*exp(j*( kix*x2+kiy*y-kfz*(z2-ac_thick(Nlay))));
        Hy_xz2=Hy_xz2+Vy*exp(j*( kix*x2+kiy*y-kfz*(z2-ac_thick(Nlay))));
        Hx_xz2=Hx_xz2+Vx*exp(j*( kix*x2+kiy*y-kfz*(z2-ac_thick(Nlay))));
        Hz_xz2=Hz_xz2+Vz*exp(j*( kix*x2+kiy*y-kfz*(z2-ac_thick(Nlay))));

        
        
% Grating region
% x-z plane
% pfEx,pfEy,pfEz,pfHx,pfHy,pfHz,pevalue
% mfEx,mfEy,mfEz,mfHx,mfHy,mfHz,mevalue

for laynt=1:Nlay

    laynt
    
pevalue=Peigvalue(1,:,laynt);
mevalue=Meigvalue(1,:,laynt);

pcnt=length(pevalue); % number of positive modes
mcnt=length(mevalue); % number of negative modes

zm=ac_thick(laynt)-lay_thick(laynt);
zp=ac_thick(laynt);
       
    zzp=zz(zzGind(laynt,1):zzGind(laynt,2))-zm;
    zzm=zz(zzGind(laynt,1):zzGind(laynt,2))-zp;
    
  for zz_cnt=1:length(zzp)
   
    for k=1:2*L
        X(k,k)   =exp(pevalue(k)*zzp(zz_cnt));
        Xinv(k,k)=exp(mevalue(k)*zzm(zz_cnt));
    end;
    
    tmp1=X*C_p(:,laynt);
    pfEx=Pf_Ex(:,:,laynt)*tmp1;
    pfEy=Pf_Ey(:,:,laynt)*tmp1;
    pfEz=Pf_Ez(:,:,laynt)*tmp1;
    pfHx=Pf_Hx(:,:,laynt)*tmp1;
    pfHy=Pf_Hy(:,:,laynt)*tmp1;
    pfHz=Pf_Hz(:,:,laynt)*tmp1;
    
    tmp2=Xinv*C_n(:,laynt);
    mfEx=Mf_Ex(:,:,laynt)*tmp2;
    mfEy=Mf_Ey(:,:,laynt)*tmp2;
    mfEz=Mf_Ez(:,:,laynt)*tmp2;
    mfHx=Mf_Hx(:,:,laynt)*tmp2;
    mfHy=Mf_Hy(:,:,laynt)*tmp2;
    mfHz=Mf_Hz(:,:,laynt)*tmp2;
    
    
     for k=1:NBx
      for l=1:NBy
          
          GEx_xz(:,zz_cnt+fc_ind(laynt,1)-1)=GEx_xz(:,zz_cnt+fc_ind(laynt,1)-1)+( pfEx((k-1)*NBy+l) + mfEx((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Ex_xz    
          GEy_xz(:,zz_cnt+fc_ind(laynt,1)-1)=GEy_xz(:,zz_cnt+fc_ind(laynt,1)-1)+( pfEy((k-1)*NBy+l) + mfEy((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Ey_xz
          GEz_xz(:,zz_cnt+fc_ind(laynt,1)-1)=GEz_xz(:,zz_cnt+fc_ind(laynt,1)-1)+( pfEz((k-1)*NBy+l) + mfEz((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Ex_xz    
          
          GHx_xz(:,zz_cnt+fc_ind(laynt,1)-1)=GHx_xz(:,zz_cnt+fc_ind(laynt,1)-1)+( pfHx((k-1)*NBy+l) + mfHx((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Hx_xz
          GHy_xz(:,zz_cnt+fc_ind(laynt,1)-1)=GHy_xz(:,zz_cnt+fc_ind(laynt,1)-1)+( pfHy((k-1)*NBy+l) + mfHy((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Hy_xz
          GHz_xz(:,zz_cnt+fc_ind(laynt,1)-1)=GHz_xz(:,zz_cnt+fc_ind(laynt,1)-1)+( pfHz((k-1)*NBy+l) + mfHz((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Hx_xz

      end; % for l
  end; % for k
      
  end; % zz_cnt


end; % for laynt
        
Ex_xz=[Ex_xz1 GEx_xz Ex_xz2];
Ey_xz=[Ey_xz1 GEy_xz Ey_xz2];
Ez_xz=[Ez_xz1 GEz_xz Ez_xz2];

Hx_xz=[Hx_xz1 GHx_xz Hx_xz2];
Hy_xz=[Hy_xz1 GHy_xz Hy_xz2];
Hz_xz=[Hz_xz1 GHz_xz Hz_xz2];

figure(1); 
subplot(3,2,1); imagesc(zz/um,xx/um,real(Ex_xz)); set(gca,'fontsize',16); colorbar; title('Ex-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); set(gca,'fontsize',16);
subplot(3,2,2); imagesc(zz/um,xx/um,real(Ey_xz)); set(gca,'fontsize',16); colorbar; title('Ey-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); set(gca,'fontsize',16);
subplot(3,2,3); imagesc(zz/um,xx/um,real(Ez_xz)); set(gca,'fontsize',16); colorbar; title('Ez-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); set(gca,'fontsize',16);
subplot(3,2,4); imagesc(zz/um,xx/um,real(Hx_xz)); set(gca,'fontsize',16); colorbar; title('Hx-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); set(gca,'fontsize',16);
subplot(3,2,5); imagesc(zz/um,xx/um,real(Hy_xz)); set(gca,'fontsize',16); colorbar; title('Hy-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); set(gca,'fontsize',16);
subplot(3,2,6); imagesc(zz/um,xx/um,real(Hz_xz)); set(gca,'fontsize',16); colorbar; title('Hz-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); set(gca,'fontsize',16);

