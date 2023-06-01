%% RCWA field visualization y-z plane
yy=5*[-Ty/2:Ty*0.002:Ty/2]';
z_start=-5*um;
z_end=5*um;
z_inc=20*nano;
zz=(z_start:z_inc:z_end);
zn=length(zz);

zrev=(z_end:-z_inc:z_start);


% y-z plane
Ex_yz=zeros(length(yy),length(zz));
Ey_yz=zeros(length(yy),length(zz));
Ez_yz=zeros(length(yy),length(zz));
Hx_yz=zeros(length(yy),length(zz));
Hy_yz=zeros(length(yy),length(zz));
Hz_yz=zeros(length(yy),length(zz));


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
 % y-z plane
  Ey_yz1=zeros(length(yy),rg1cnt);
  Ex_yz1=zeros(length(yy),rg1cnt);
  Ez_yz1=zeros(length(yy),rg1cnt);
  Hy_yz1=zeros(length(yy),rg1cnt);
  Hx_yz1=zeros(length(yy),rg1cnt);
  Hz_yz1=zeros(length(yy),rg1cnt);     

 % Grating region
 
  GEy_yz=zeros(length(yy),sum(rgGcnt));
  GEx_yz=zeros(length(yy),sum(rgGcnt));
  GEz_yz=zeros(length(yy),sum(rgGcnt));
  GHy_yz=zeros(length(yy),sum(rgGcnt));
  GHx_yz=zeros(length(yy),sum(rgGcnt));
  GHz_yz=zeros(length(yy),sum(rgGcnt));     

  % region 2
  
  Ey_yz2=zeros(length(yy),rg2cnt);
  Ex_yz2=zeros(length(yy),rg2cnt);
  Ez_yz2=zeros(length(yy),rg2cnt);
  Hy_yz2=zeros(length(yy),rg2cnt);
  Hx_yz2=zeros(length(yy),rg2cnt);
  Hz_yz2=zeros(length(yy),rg2cnt);     

  Iden=zeros(2*L,2*L);
for k=1:2*L
Iden(k,k)=1;   
end;
          
       
%% .    input field - plane wave

input_pol_mode='TM';
%input_pol_mode='TE';

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
Hinc=-Vh*Einc;
Vy=Hinc(centx);
Vx=Hinc(L+centx);
Vz=(kx0*Vx+ky0*Vy)/kz0;

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
Hinc=-Vh*Einc;
Vy=Hinc(centx);
Vx=Hinc(L+centx);
Vz=(kx0*Vx+ky0*Vy)/kz0;

te_Vy=Vy;
te_Vx=Vx;
te_Vz=Vz;

end;

ET1_temp=TTb*Einc;
HT1_temp=-Vh*ET1_temp;
ET2_temp=RRb*Einc;
HT2_temp=Vh*ET2_temp;


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
      
      ET1y_cof(k,l)=ET1_y((k-1)*NBy+l );
      ET1x_cof(k,l)=ET1_x((k-1)*NBy+l );
	  ET1z_cof(k,l)=( kx_vc(k)*ET1x_cof(k,l)+ky_vc(l)*ET1y_cof(k,l) )/kz_vc(k,l); 
      HT1y_cof(k,l)=HT1_y((k-1)*NBy+l );
      HT1x_cof(k,l)=HT1_x((k-1)*NBy+l );
	  HT1z_cof(k,l)=( kx_vc(k)*HT1x_cof(k,l)+ky_vc(l)*HT1y_cof(k,l) )/kz_vc(k,l); 
       
      ET2y_cof(k,l)=ET2_y((k-1)*NBy+l);
      ET2x_cof(k,l)=ET2_x((k-1)*NBy+l );
      ET2z_cof(k,l)=-( kx_vc(k)*ET2x_cof(k,l)+ky_vc(l)*ET2y_cof(k,l) )/kz_vc(k,l); 
      HT2y_cof(k,l)=HT2_y((k-1)*NBy+l);
      HT2x_cof(k,l)=HT2_x((k-1)*NBy+l );
      HT2z_cof(k,l)=-( kx_vc(k)*HT2x_cof(k,l)+ky_vc(l)*HT2y_cof(k,l) )/kz_vc(k,l); 
   
   end;
end;



% diffraction efficiency

DEt1=zeros(NBx,NBy);
DEt2=zeros(NBx,NBy);
for k=1:NBx
   for l=1:NBy
      
      DEt1(k,l)=abs(abs(ET1x_cof(k,l))^2+abs(ET1y_cof(k,l))^2+abs(ET1z_cof(k,l))^2)*real(kz_vc(k,l)/(kz0));
      DEt2(k,l)=abs(abs(ET2x_cof(k,l))^2+abs(ET2y_cof(k,l))^2+abs(ET2z_cof(k,l))^2)*real(kz_vc(k,l)/(kz0));
      
   end;   
end;

total_energy=sum(sum(DEt1))+sum(sum(DEt2));

% y-z Field visualization
x=0;
% region 1
[z1 y1]=meshgrid(zz(zz1ind(1):zz1ind(2)),yy); % region 1
  for k=1:NBx
           for l=1:NBy
           Ey_yz1=Ey_yz1+ET1y_cof(k,l)*exp(j*( kx_vc(k)*x+ky_vc(l)*y1-kz_vc(k,l)*z1 ));   
           Ex_yz1=Ex_yz1+ET1x_cof(k,l)*exp(j*( kx_vc(k)*x+ky_vc(l)*y1-kz_vc(k,l)*z1 ));
           Ez_yz1=Ez_yz1+ET1z_cof(k,l)*exp(j*( kx_vc(k)*x+ky_vc(l)*y1-kz_vc(k,l)*z1 ));
           
           Hy_yz1=Hy_yz1+HT1y_cof(k,l)*exp(j*( kx_vc(k)*x+ky_vc(l)*y1-kz_vc(k,l)*z1 ));   
           Hx_yz1=Hx_yz1+HT1x_cof(k,l)*exp(j*( kx_vc(k)*x+ky_vc(l)*y1-kz_vc(k,l)*z1 ));
           Hz_yz1=Hz_yz1+HT1z_cof(k,l)*exp(j*( kx_vc(k)*x+ky_vc(l)*y1-kz_vc(k,l)*z1 ));
           end;
        end;
        
%        Ey_yz1=Ey_yz1+Uy*exp(j*( kx0*x+ky0*y1+kz0*z1));
%        Ex_yz1=Ex_yz1+Ux*exp(j*( kx0*x+ky0*y1+kz0*z1));
%        Ez_yz1=Ez_yz1+Uz*exp(j*( kx0*x+ky0*y1+kz0*z1));
%        Hy_yz1=Hy_yz1+Vy*exp(j*( kx0*x+ky0*y1+kz0*z1));
%        Hx_yz1=Hx_yz1+Vx*exp(j*( kx0*x+ky0*y1+kz0*z1));
%        Hz_yz1=Hz_yz1+Vz*exp(j*( kx0*x+ky0*y1+kz0*z1));


% region 3

[z2 y2]=meshgrid(zz(zz2ind(1):zz2ind(2)),yy); % region 2
    for k=1:NBx
         for l=1:NBy
           Ey_yz2=Ey_yz2+ET2y_cof(k,l)*exp(j*kz_vc(k,l)*(z2-ac_thick(Nlay))).*exp(j*( kx_vc(k)*x+ky_vc(l)*y2 ));   
           Ex_yz2=Ex_yz2+ET2x_cof(k,l)*exp(j*kz_vc(k,l)*(z2-ac_thick(Nlay))).*exp(j*( kx_vc(k)*x+ky_vc(l)*y2 ));
           Ez_yz2=Ez_yz2+ET2z_cof(k,l)*exp(j*kz_vc(k,l)*(z2-ac_thick(Nlay))).*exp(j*( kx_vc(k)*x+ky_vc(l)*y2 ));
           Hy_yz2=Hy_yz2+HT2y_cof(k,l)*exp(j*kz_vc(k,l)*(z2-ac_thick(Nlay))).*exp(j*( kx_vc(k)*x+ky_vc(l)*y2 ));   
           Hx_yz2=Hx_yz2+HT2x_cof(k,l)*exp(j*kz_vc(k,l)*(z2-ac_thick(Nlay))).*exp(j*( kx_vc(k)*x+ky_vc(l)*y2 ));
           Hz_yz2=Hz_yz2+HT2z_cof(k,l)*exp(j*kz_vc(k,l)*(z2-ac_thick(Nlay))).*exp(j*( kx_vc(k)*x+ky_vc(l)*y2 ));

            end;
        end;
        
        Ey_yz2=Ey_yz2+Uy*exp(j*( kx0*x+ky0*y2-kz0*(z2-ac_thick(Nlay))));
        Ex_yz2=Ex_yz2+Ux*exp(j*( kx0*x+ky0*y2-kz0*(z2-ac_thick(Nlay))));
        Ez_yz2=Ez_yz2+Uz*exp(j*( kx0*x+ky0*y2-kz0*(z2-ac_thick(Nlay))));
        Hy_yz2=Hy_yz2+Vy*exp(j*( kx0*x+ky0*y2-kz0*(z2-ac_thick(Nlay))));
        Hx_yz2=Hx_yz2+Vx*exp(j*( kx0*x+ky0*y2-kz0*(z2-ac_thick(Nlay))));
        Hz_yz2=Hz_yz2+Vz*exp(j*( kx0*x+ky0*y2-kz0*(z2-ac_thick(Nlay))));
        
        %         
% region 2 : Grating region
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
          
          GEx_yz(:,zz_cnt+fc_ind(laynt,1)-1)=GEx_yz(:,zz_cnt+fc_ind(laynt,1)-1)+( pfEx((k-1)*NBy+l) + mfEx((k-1)*NBy+l) )*exp(j*( kx_vc(k)*x+ky_vc(l)*yy));   % Ex_xz    
          GEy_yz(:,zz_cnt+fc_ind(laynt,1)-1)=GEy_yz(:,zz_cnt+fc_ind(laynt,1)-1)+( pfEy((k-1)*NBy+l) + mfEy((k-1)*NBy+l) )*exp(j*( kx_vc(k)*x+ky_vc(l)*yy));   % Ey_xz
          GEz_yz(:,zz_cnt+fc_ind(laynt,1)-1)=GEz_yz(:,zz_cnt+fc_ind(laynt,1)-1)+( pfEz((k-1)*NBy+l) + mfEz((k-1)*NBy+l) )*exp(j*( kx_vc(k)*x+ky_vc(l)*yy));   % Ex_xz    
          
          GHx_yz(:,zz_cnt+fc_ind(laynt,1)-1)=GHx_yz(:,zz_cnt+fc_ind(laynt,1)-1)+( pfHx((k-1)*NBy+l) + mfHx((k-1)*NBy+l) )*exp(j*( kx_vc(k)*x+ky_vc(l)*yy));   % Hx_xz
          GHy_yz(:,zz_cnt+fc_ind(laynt,1)-1)=GHy_yz(:,zz_cnt+fc_ind(laynt,1)-1)+( pfHy((k-1)*NBy+l) + mfHy((k-1)*NBy+l) )*exp(j*( kx_vc(k)*x+ky_vc(l)*yy));   % Hy_xz
          GHz_yz(:,zz_cnt+fc_ind(laynt,1)-1)=GHz_yz(:,zz_cnt+fc_ind(laynt,1)-1)+( pfHz((k-1)*NBy+l) + mfHz((k-1)*NBy+l) )*exp(j*( kx_vc(k)*x+ky_vc(l)*yy));   % Hx_xz

      end; % for l
  end; % for k
      
  end; % zz_cnt

end; % for laynt   

Ex_yz=[Ex_yz1 GEx_yz  Ex_yz2];
Ey_yz=[Ey_yz1 GEy_yz  Ey_yz2];
Ez_yz=[Ez_yz1 GEz_yz  Ez_yz2];

Hx_yz=[Hx_yz1 GHx_yz  Hx_yz2];
Hy_yz=[Hy_yz1 GHy_yz  Hy_yz2];
Hz_yz=[Hz_yz1 GHz_yz  Hz_yz2];

figure(1); 
subplot(3,2,1); imagesc(zz/um,yy/um,real(Ex_yz)); set(gca,'fontsize',16); colorbar; title('Ex-yz'); xlabel('z-axis(um)');ylabel('y-axis(um)'); set(gca,'fontsize',16);
subplot(3,2,2); imagesc(zz/um,yy/um,real(Ey_yz)); set(gca,'fontsize',16); colorbar; title('Ey-yz'); xlabel('z-axis(um)');ylabel('y-axis(um)'); set(gca,'fontsize',16);
subplot(3,2,3); imagesc(zz/um,yy/um,real(Ez_yz)); set(gca,'fontsize',16); colorbar; title('Ez-yz'); xlabel('z-axis(um)');ylabel('y-axis(um)'); set(gca,'fontsize',16);
subplot(3,2,4); imagesc(zz/um,yy/um,real(Hx_yz)); set(gca,'fontsize',16); colorbar; title('Hx-yz'); xlabel('z-axis(um)');ylabel('y-axis(um)'); set(gca,'fontsize',16);
subplot(3,2,5); imagesc(zz/um,yy/um,real(Hy_yz)); set(gca,'fontsize',16); colorbar; title('Hy-yz'); xlabel('z-axis(um)');ylabel('y-axis(um)'); set(gca,'fontsize',16);
subplot(3,2,6); imagesc(zz/um,yy/um,real(Hz_yz)); set(gca,'fontsize',16); colorbar; title('Hz-yz'); xlabel('z-axis(um)');ylabel('y-axis(um)'); set(gca,'fontsize',16);


