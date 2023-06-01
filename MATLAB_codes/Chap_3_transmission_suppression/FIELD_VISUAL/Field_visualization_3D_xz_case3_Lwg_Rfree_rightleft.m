%% RCWA field visualization x-z plane
% right to left

xx=5*[-Tx/2:Tx*0.002:Tx/2]';
z_start=-10*um;
z_end=20*um;
z_inc=50*nano;
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
rg2cnt=zeros(Nlay,1);   % region 2 : grating region 
rg3cnt=0;               % region 3

for znt=1:zn
    
    zv=zz(znt);
    
    if zv <= 0
      rg1cnt=rg1cnt+1;
    end;
    
    if (0 < zv) & (zv <= ac_thick(Nlay))
            
       for laynt=1:Nlay
             if  ( ( ac_thick(laynt)-lay_thick(laynt) ) < zv ) & ( zv <= ac_thick(laynt)  )
              rg2cnt(laynt)=rg2cnt(laynt)+1; 
             end;
       end;
    
    end;
    
     if ( ac_thick(Nlay) < zv )
        rg3cnt=rg3cnt+1;
    end;
    
    
end; % for znt

 % region 1
 zz1ind(1,1)=1;
 zz1ind(1,2)=rg1cnt;
 
 % region 2 : grating region
 for laynt=1:Nlay
     if laynt==1
 zz2ind(laynt,1)=rg1cnt+1;
 zz2ind(laynt,2)=zz2ind(laynt,1)+rg2cnt(laynt)-1 ;
 
     else
 zz2ind(laynt,1)=zz2ind(laynt-1,2)+1;
 zz2ind(laynt,2)=zz2ind(laynt,1)+rg2cnt(laynt)-1 ;
     end;
 end;
 
 fc_ind=zz2ind-rg1cnt; % field calculation index for region 2
 
 % region 3
 zz3ind(1,1)=rg1cnt+sum(rg2cnt)+1;
 zz3ind(1,2)=rg1cnt+sum(rg2cnt)+rg3cnt;
 
 
 % region 1
 % x-z plane
  Ey_xz1=zeros(length(xx),rg1cnt);
  Ex_xz1=zeros(length(xx),rg1cnt);
  Ez_xz1=zeros(length(xx),rg1cnt);
  Hy_xz1=zeros(length(xx),rg1cnt);
  Hx_xz1=zeros(length(xx),rg1cnt);
  Hz_xz1=zeros(length(xx),rg1cnt);     

 % region 2
 
  Ey_xz2=zeros(length(xx),sum(rg2cnt));
  Ex_xz2=zeros(length(xx),sum(rg2cnt));
  Ez_xz2=zeros(length(xx),sum(rg2cnt));
  Hy_xz2=zeros(length(xx),sum(rg2cnt));
  Hx_xz2=zeros(length(xx),sum(rg2cnt));
  Hz_xz2=zeros(length(xx),sum(rg2cnt));     

  % region 3
  
  Ey_xz3=zeros(length(xx),rg3cnt);
  Ex_xz3=zeros(length(xx),rg3cnt);
  Ez_xz3=zeros(length(xx),rg3cnt);
  Hy_xz3=zeros(length(xx),rg3cnt);
  Hx_xz3=zeros(length(xx),rg3cnt);
  Hz_xz3=zeros(length(xx),rg3cnt);     

  Iden=zeros(2*L,2*L);
for k=1:2*L
Iden(k,k)=1;   
end;
          
       
%% .    input field - plane wave

% Uy=tm_Uy;
% Ux=tm_Ux;
% Uz=tm_Uz;
% 
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

% tm_Vy=Vy;
% tm_Vx=Vx;
% tm_Vz=Vz;

te_Vy=Vy;
te_Vx=Vx;
te_Vz=Vz;

% region I
CT1=TTb*Einc;
Lay_cof1=zeros(4*L,1);
Lay_cof1(1:2*L,1)=0;
Lay_cof1(2*L+1:4*L,1)=CT1;


% region III
ET3_temp=RRb*Einc;
HT3_temp=VhII*ET3_temp;

% grating
for laynt=1:Nlay
C_temp(:,laynt)=Cb(:,:,laynt)*Einc;
C_p(:,laynt)=C_temp(1:2*L,laynt);          % Positive coupling coefficients
C_n(:,laynt)=C_temp(2*L+1:4*L,laynt);      % Negative coupling coefficients  
end;

ET3_y=zeros(L,1);
ET3_x=zeros(L,1);
HT3_y=zeros(L,1);
HT3_x=zeros(L,1);

ET3_y=ET3_temp(1:L);
ET3_x=ET3_temp(L+1:2*L);
HT3_y=HT3_temp(1:L);
HT3_x=HT3_temp(L+1:2*L);

ET3x_cof=zeros(NBx,NBy);
ET3y_cof=zeros(NBx,NBy);
ET3z_cof=zeros(NBx,NBy);
HT3x_cof=zeros(NBx,NBy);
HT3y_cof=zeros(NBx,NBy);
HT3z_cof=zeros(NBx,NBy);

for k=1:NBx
   for l=1:NBy
      
      ET3y_cof(k,l)=ET3_y((k-1)*NBy+l);
      ET3x_cof(k,l)=ET3_x((k-1)*NBy+l );
      ET3z_cof(k,l)=-( kx_tra(k)*ET3x_cof(k,l)+ky_tra(l)*ET3y_cof(k,l) )/kz_tra(k,l); 
      HT3y_cof(k,l)=HT3_y((k-1)*NBy+l);
      HT3x_cof(k,l)=HT3_x((k-1)*NBy+l );
      HT3z_cof(k,l)=-( kx_tra(k)*HT3x_cof(k,l)+ky_tra(l)*HT3y_cof(k,l) )/kz_tra(k,l); 

   
   end;
end;

% 
% % diffraction efficiency
% 
% DEt1=zeros(NBx,NBy);
% DEt3=zeros(NBx,NBy);
% for k=1:NBx
%    for l=1:NBy
%       
%       DEt1(k,l)=abs(abs(ET1x_cof(k,l))^2+abs(ET1y_cof(k,l))^2+abs(ET1z_cof(k,l))^2)*real(kz_ref(k,l)/(kiz));
%       DEt3(k,l)=abs(abs(ET3x_cof(k,l))^2+abs(ET3y_cof(k,l))^2+abs(ET3z_cof(k,l))^2)*real(kz_tra(k,l)/(kiz));
%       
%    end;   
% end;
%    
% 
% total_energy=sum(sum(DEt1))+sum(sum(DEt3));

% x-z Field visualization
y=0;
% region 1

laynt=1;
    
pevalue=Peigvalue(1,:,laynt);
mevalue=Meigvalue(1,:,laynt);

pcnt=length(pevalue); % number of positive modes
mcnt=length(mevalue); % number of negative modes

zm=ac_thick(laynt)-lay_thick(laynt);
zp=ac_thick(laynt);
       
    zzp=zz(zz1ind(1):zz1ind(2))-zm;
    zzm=zz(zz1ind(1):zz1ind(2))-zp;
    
  for zz_cnt=1:length(zzp)
   
    for k=1:2*L
        X(k,k)   =exp(pevalue(k)*zzp(zz_cnt));
        Xinv(k,k)=exp(mevalue(k)*zzm(zz_cnt)); % there is no reflection mode
    end;
    
    tmp1=X*Lay_cof1(1:2*L,1);
    pfEx=Pf_Ex(:,:,laynt)*tmp1;
    pfEy=Pf_Ey(:,:,laynt)*tmp1;
    pfEz=Pf_Ez(:,:,laynt)*tmp1;
    pfHx=Pf_Hx(:,:,laynt)*tmp1;
    pfHy=Pf_Hy(:,:,laynt)*tmp1;
    pfHz=Pf_Hz(:,:,laynt)*tmp1;
    
    tmp2=Xinv*Lay_cof1(2*L+1:4*L,1);
    mfEx=Mf_Ex(:,:,laynt)*tmp2;
    mfEy=Mf_Ey(:,:,laynt)*tmp2;
    mfEz=Mf_Ez(:,:,laynt)*tmp2;
    mfHx=Mf_Hx(:,:,laynt)*tmp2;
    mfHy=Mf_Hy(:,:,laynt)*tmp2;
    mfHz=Mf_Hz(:,:,laynt)*tmp2;
    
    
     for k=1:NBx
      for l=1:NBy
          
          Ex_xz1(:,zz_cnt)=Ex_xz1(:,zz_cnt)+( pfEx((k-1)*NBy+l) + mfEx((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Ex_xz    
          Ey_xz1(:,zz_cnt)=Ey_xz1(:,zz_cnt)+( pfEy((k-1)*NBy+l) + mfEy((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Ey_xz
          Ez_xz1(:,zz_cnt)=Ez_xz1(:,zz_cnt)+( pfEz((k-1)*NBy+l) + mfEz((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Ex_xz    
          
          Hx_xz1(:,zz_cnt)=Hx_xz1(:,zz_cnt)+( pfHx((k-1)*NBy+l) + mfHx((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Hx_xz
          Hy_xz1(:,zz_cnt)=Hy_xz1(:,zz_cnt)+( pfHy((k-1)*NBy+l) + mfHy((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Hy_xz
          Hz_xz1(:,zz_cnt)=Hz_xz1(:,zz_cnt)+( pfHz((k-1)*NBy+l) + mfHz((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Hx_xz

      end; % for l
  end; % for k
      
  end; % zz_cnt
  
  
% region 3

[z3 x3]=meshgrid(zz(zz3ind(1):zz3ind(2)),xx); % region 3
    for k=1:NBx
         for l=1:NBy
           Ey_xz3=Ey_xz3+ET3y_cof(k,l)*exp(j*kz_tra(k,l)*(z3-ac_thick(Nlay))).*exp(j*( kx_tra(k)*x3+ky_tra(l)*y ));   
           Ex_xz3=Ex_xz3+ET3x_cof(k,l)*exp(j*kz_tra(k,l)*(z3-ac_thick(Nlay))).*exp(j*( kx_tra(k)*x3+ky_tra(l)*y ));
           Ez_xz3=Ez_xz3+ET3z_cof(k,l)*exp(j*kz_tra(k,l)*(z3-ac_thick(Nlay))).*exp(j*( kx_tra(k)*x3+ky_tra(l)*y ));
           Hy_xz3=Hy_xz3+HT3y_cof(k,l)*exp(j*kz_tra(k,l)*(z3-ac_thick(Nlay))).*exp(j*( kx_tra(k)*x3+ky_tra(l)*y ));   
           Hx_xz3=Hx_xz3+HT3x_cof(k,l)*exp(j*kz_tra(k,l)*(z3-ac_thick(Nlay))).*exp(j*( kx_tra(k)*x3+ky_tra(l)*y ));
           Hz_xz3=Hz_xz3+HT3z_cof(k,l)*exp(j*kz_tra(k,l)*(z3-ac_thick(Nlay))).*exp(j*( kx_tra(k)*x3+ky_tra(l)*y ));

          end;
     end;
     
        Ey_xz3=Ey_xz3+Uy*exp(j*( kix*x3+kiy*y-kfz*(z3-ac_thick(Nlay))));
        Ex_xz3=Ex_xz3+Ux*exp(j*( kix*x3+kiy*y-kfz*(z3-ac_thick(Nlay))));
        Ez_xz3=Ez_xz3+Uz*exp(j*( kix*x3+kiy*y-kz0*(z3-ac_thick(Nlay))));
        Hy_xz3=Hy_xz3+Vy*exp(j*( kix*x3+kiy*y-kz0*(z3-ac_thick(Nlay))));
        Hx_xz3=Hx_xz3+Vx*exp(j*( kix*x3+kiy*y-kz0*(z3-ac_thick(Nlay))));
        Hz_xz3=Hz_xz3+Vz*exp(j*( kix*x3+kiy*y-kz0*(z3-ac_thick(Nlay))));

 
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
       
    zzp=zz(zz2ind(laynt,1):zz2ind(laynt,2))-zm;
    zzm=zz(zz2ind(laynt,1):zz2ind(laynt,2))-zp;
    
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
          
          Ex_xz2(:,zz_cnt+fc_ind(laynt,1)-1)=Ex_xz2(:,zz_cnt+fc_ind(laynt,1)-1)+( pfEx((k-1)*NBy+l) + mfEx((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Ex_xz    
          Ey_xz2(:,zz_cnt+fc_ind(laynt,1)-1)=Ey_xz2(:,zz_cnt+fc_ind(laynt,1)-1)+( pfEy((k-1)*NBy+l) + mfEy((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Ey_xz
          Ez_xz2(:,zz_cnt+fc_ind(laynt,1)-1)=Ez_xz2(:,zz_cnt+fc_ind(laynt,1)-1)+( pfEz((k-1)*NBy+l) + mfEz((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Ex_xz    
          
          Hx_xz2(:,zz_cnt+fc_ind(laynt,1)-1)=Hx_xz2(:,zz_cnt+fc_ind(laynt,1)-1)+( pfHx((k-1)*NBy+l) + mfHx((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Hx_xz
          Hy_xz2(:,zz_cnt+fc_ind(laynt,1)-1)=Hy_xz2(:,zz_cnt+fc_ind(laynt,1)-1)+( pfHy((k-1)*NBy+l) + mfHy((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Hy_xz
          Hz_xz2(:,zz_cnt+fc_ind(laynt,1)-1)=Hz_xz2(:,zz_cnt+fc_ind(laynt,1)-1)+( pfHz((k-1)*NBy+l) + mfHz((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Hx_xz

      end; % for l
  end; % for k
      
  end; % zz_cnt


end; % for laynt
        
Ex_xz=[Ex_xz1 Ex_xz2 Ex_xz3];
Ey_xz=[Ey_xz1 Ey_xz2 Ey_xz3];
Ez_xz=[Ez_xz1 Ez_xz2 Ez_xz3];

Hx_xz=[Hx_xz1 Hx_xz2 Hx_xz3];
Hy_xz=[Hy_xz1 Hy_xz2 Hy_xz3];
Hz_xz=[Hz_xz1 Hz_xz2 Hz_xz3];


figure(1); imagesc(zz/um,xx/um,real(Ex_xz)); set(gca,'fontsize',16); colorbar; title('Ex-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)');colormap(hot); set(gca,'fontsize',16);
figure(2); imagesc(zz/um,xx/um,real(Ey_xz)); set(gca,'fontsize',16); colorbar; title('Ey-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)');colormap(hot); set(gca,'fontsize',16);
figure(3); imagesc(zz/um,xx/um,real(Ez_xz)); set(gca,'fontsize',16); colorbar; title('Ez-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)');colormap(hot); set(gca,'fontsize',16);


figure(4); imagesc(zz/um,xx/um,real(Hx_xz)); set(gca,'fontsize',16); colorbar; title('Hx-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);
figure(5); imagesc(zz/um,xx/um,real(Hy_xz)); set(gca,'fontsize',16); colorbar; title('Hy-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);
figure(6); imagesc(zz/um,xx/um,real(Hz_xz)); set(gca,'fontsize',16); colorbar; title('Hz-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);


% %% save results ( fig, png, mat for Ex_xz Ey_xz Ez_xz )
% %  save results (fig, png, mat, for Hx_xz Hy_xz Hz_xz)
% close all;
% 
% width_thick = sprintf('_thick_%04d_width_%04d',film_thick_nm,Wx_nm);
% fname_Ex_xz  = sprintf('data/Ex/%s%s','Ex_xz',width_thick);
% fname_Ey_xz  = sprintf('data/Ey/%s%s','Ey_xz',width_thick);
% fname_Ez_xz  = sprintf('data/Ez/%s%s','Ez_xz',width_thick);
% fname_Hx_xz  = sprintf('data/Hx/%s%s','Hx_xz',width_thick);
% fname_Hy_xz  = sprintf('data/Hy/%s%s','Hy_xz',width_thick);
% fname_Hz_xz  = sprintf('data/Hz/%s%s','Hz_xz',width_thick);
% 
% fname_Ex_xz_p  = sprintf('data/Ex/%s%s.png','Ex_xz',width_thick);  % for png
% fname_Ey_xz_p  = sprintf('data/Ey/%s%s.png','Ey_xz',width_thick);
% fname_Ez_xz_p  = sprintf('data/Ez/%s%s.png','Ez_xz',width_thick);
% fname_Hx_xz_p  = sprintf('data/Hx/%s%s.png','Hx_xz',width_thick);
% fname_Hy_xz_p  = sprintf('data/Hy/%s%s.png','Hy_xz',width_thick);
% fname_Hz_xz_p  = sprintf('data/Hz/%s%s.png','Hz_xz',width_thick);
% 
% % Ex_xz
% save(fname_Ex_xz,'Ex_xz');
% figure(1);
% set(gca,'fontsize',16);
% imagesc(xx/um,zz/um,abs(Ex_xz)');
% caxis([0 2]);
% axis equal;
% axis([xx(1)/um xx(end)/um zz(1)/um zz(end)/um]);
% set(gca,'ydir','normal');
% xlabel('x [\mum]');
% ylabel('z [\mum]');
% title(sprintf('|E_x| (xz-plane)     t : %04d nm    w : %04d nm',film_thick_nm,Wx_nm));
% colorbar;
% set(gca,'fontsize',16);
% %               figfs
% hgsave(fname_Ex_xz);
% %               png
% print('-dpng',fname_Ex_xz_p);
% % Ex_xz_map = imread(fname_Ex_xz_p);
% % new_Ex_xz_map=Ex_xz_map(310:580,50:1150,:);
% % imwrite(new_Ex_xz_map,fname_Ex_xz_p,'png');
% %
% 
% % Ey_xz
% save(fname_Ey_xz,'Ey_xz');
% figure(2);
% set(gca,'fontsize',16);
% imagesc(xx/um,zz/um,abs(Ey_xz)');
% caxis([0 2]);
% axis equal;
% axis([xx(1)/um xx(end)/um zz(1)/um zz(end)/um]);
% set(gca,'ydir','normal');
% xlabel('x [\mum]');
% ylabel('z [\mum]');
% title(sprintf('|E_y| (xz-plane)     t : %04d nm    w : %04d nm',film_thick_nm,Wx_nm));
% colorbar;
% set(gca,'fontsize',16);
% %               figfs
% hgsave(fname_Ey_xz);
% %               png
% print('-dpng',fname_Ey_xz_p);
% % Ey_xz_map = imread(fname_Ey_xz_p);
% % new_Ey_xz_map=Ey_xz_map(310:580,50:1150,:);
% % imwrite(new_Ey_xz_map,fname_Ey_xz_p,'png');
% %
% 
% % Ez_xz
% save(fname_Ez_xz,'Ez_xz');
% figure(3);
% set(gca,'fontsize',16);
% imagesc(xx/um,zz/um,abs(Ez_xz)');
% caxis([0 2]);
% axis equal;
% axis([xx(1)/um xx(end)/um zz(1)/um zz(end)/um]);
% set(gca,'ydir','normal');
% xlabel('x [\mum]');
% ylabel('z [\mum]');
% title(sprintf('|E_z| (xz-plane)     t : %04d nm    w : %04d nm',film_thick_nm,Wx_nm));
% colorbar;
% set(gca,'fontsize',16);
% %               figfs
% hgsave(fname_Ez_xz);
% %               png
% print('-dpng',fname_Ez_xz_p);
% % Ez_xz_map = imread(fname_Ez_xz_p);
% % new_Ez_xz_map=Ez_xz_map(310:580,50:1150,:);
% % imwrite(new_Ez_xz_map,fname_Ez_xz_p,'png');
% %
% 
% % Hx_xz
% save(fname_Hx_xz,'Hx_xz');
% figure(4);
% set(gca,'fontsize',16);
% imagesc(xx/um,zz/um,abs(Hx_xz)');
% caxis([0 2*1e-2]);
% axis equal;
% axis([xx(1)/um xx(end)/um zz(1)/um zz(end)/um]);
% set(gca,'ydir','normal');
% xlabel('x [\mum]');
% ylabel('z [\mum]');
% title(sprintf('|H_x| (xz-plane)     t : %04d nm    w : %04d nm',film_thick_nm,Wx_nm));
% colorbar;
% set(gca,'fontsize',16);
% %               figfs
% hgsave(fname_Hx_xz);
% %               png
% print('-dpng',fname_Hx_xz_p);
% % Hx_xz_map = imread(fname_Hx_xz_p);
% % new_Hx_xz_map=Hx_xz_map(310:580,50:1150,:);
% % imwrite(new_Hx_xz_map,fname_Hx_xz_p,'png');
% % %
% 
% % Hy_xz
% save(fname_Hy_xz,'Hy_xz');
% figure(5);
% set(gca,'fontsize',16);
% imagesc(xx/um,zz/um,abs(Hy_xz)');
% caxis([0 2*1e-2]);
% axis equal;
% axis([xx(1)/um xx(end)/um zz(1)/um zz(end)/um]);
% set(gca,'ydir','normal');
% xlabel('x [\mum]');
% ylabel('z [\mum]');
% title(sprintf('|H_y| (xz-plane)     t : %04d nm    w : %04d nm',film_thick_nm,Wx_nm));
% colorbar;
% set(gca,'fontsize',16);
% %               figfs
% hgsave(fname_Hy_xz);
% %               png
% print('-dpng',fname_Hy_xz_p);
% % Hy_xz_map = imread(fname_Hy_xz_p);
% % new_Hy_xz_map=Hy_xz_map(310:580,50:1150,:);
% % imwrite(new_Hy_xz_map,fname_Hy_xz_p,'png');
% %
% 
% % Hz_xz
% save(fname_Hz_xz,'Hz_xz');
% figure(6);
% set(gca,'fontsize',16);
% imagesc(xx/um,zz/um,abs(Hz_xz)');
% caxis([0 2*1e-2]);
% axis equal;
% axis([xx(1)/um xx(end)/um zz(1)/um zz(end)/um]);
% set(gca,'ydir','normal');
% xlabel('x [\mum]');
% ylabel('z [\mum]');
% title(sprintf('|H_z| (xz-plane)     t : %04d nm    w : %04d nm',film_thick_nm,Wx_nm));
% colorbar;
% set(gca,'fontsize',16);
% %               figfs
% hgsave(fname_Hz_xz);
% %               png
% print('-dpng',fname_Hz_xz_p);
% % Hz_xz_map = imread(fname_Hz_xz_p);
% % new_Hz_xz_map=Hz_xz_map(310:580,50:1150,:);
% % imwrite(new_Hz_xz_map,fname_Hz_xz_p,'png');
% % %

