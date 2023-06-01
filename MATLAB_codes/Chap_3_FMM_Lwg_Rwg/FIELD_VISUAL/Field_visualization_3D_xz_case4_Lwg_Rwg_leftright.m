%% Field visualization case 4 left-right
% left to right

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
       
%% .    input field - eigenmode of left waveguide

% excitation eigenmode
Einc=zeros(2*L,1);
laynt=1; % left waveguide
pevalue=Peigvalue(1,:,laynt);
mevalue=Meigvalue(1,:,laynt);
[p_ind1 p_ind2]=sort( -real(pevalue) );
[m_ind1 m_ind2]=sort(  real(mevalue) );
% positive eigenmode  
Einc(p_ind2(1))=1;

% region 1
CR1=RRa*Einc;
Lay_cof1=zeros(4*L,1);
Lay_cof1(1:2*L,1)=Einc;
Lay_cof1(2*L+1:4*L,1)=CR1;

% region 2
CT2=TTa*Einc;
Lay_cof2=zeros(4*L,1);
Lay_cof2(1:2*L,1)=CT2;

% grating
for laynt=1:Nlay
C_temp(:,laynt)=Ca(:,:,laynt)*Einc;
C_p(:,laynt)=C_temp(1:2*L,laynt);          % Positive coupling coefficients
C_n(:,laynt)=C_temp(2*L+1:4*L,laynt);      % Negative coupling coefficients  
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
        
        if abs(Xinv(k,k)) == Inf
          Xinv(k,k)=0;              %% precision error
        end;
        if abs(X(k,k)) == Inf
            X(k,k)=0;               %% precision error
        end;
        
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
  
  
% region 2
laynt=Nlay;
    
pevalue=Peigvalue(1,:,laynt);
mevalue=Meigvalue(1,:,laynt);

pcnt=length(pevalue); % number of positive modes
mcnt=length(mevalue); % number of negative modes

zm=ac_thick(laynt)-lay_thick(laynt);
zp=ac_thick(laynt);
       
    zzp=zz(zz2ind(1):zz2ind(2))-zm;
    zzm=zz(zz2ind(1):zz2ind(2))-zp;
    
  for zz_cnt=1:length(zzp)
   
    for k=1:2*L
        X(k,k)   =exp(pevalue(k)*zzp(zz_cnt));
        Xinv(k,k)=exp(mevalue(k)*zzm(zz_cnt)); % there is no reflection mode
        if abs(Xinv(k,k)) == Inf
          Xinv(k,k)=0;              %% precision error
        end;
        if abs(X(k,k)) == Inf
            X(k,k)=0;               %% precision error
        end;
    end;
    
    tmp1=X*Lay_cof2(1:2*L,1);
    pfEx=Pf_Ex(:,:,laynt)*tmp1;
    pfEy=Pf_Ey(:,:,laynt)*tmp1;
    pfEz=Pf_Ez(:,:,laynt)*tmp1;
    pfHx=Pf_Hx(:,:,laynt)*tmp1;
    pfHy=Pf_Hy(:,:,laynt)*tmp1;
    pfHz=Pf_Hz(:,:,laynt)*tmp1;
    
    tmp2=Xinv*Lay_cof2(2*L+1:4*L,1);
    mfEx=Mf_Ex(:,:,laynt)*tmp2;
    mfEy=Mf_Ey(:,:,laynt)*tmp2;
    mfEz=Mf_Ez(:,:,laynt)*tmp2;
    mfHx=Mf_Hx(:,:,laynt)*tmp2;
    mfHy=Mf_Hy(:,:,laynt)*tmp2;
    mfHz=Mf_Hz(:,:,laynt)*tmp2;
    
    
     for k=1:NBx
      for l=1:NBy
          
          Ex_xz2(:,zz_cnt)=Ex_xz2(:,zz_cnt)+( pfEx((k-1)*NBy+l) + mfEx((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Ex_xz    
          Ey_xz2(:,zz_cnt)=Ey_xz2(:,zz_cnt)+( pfEy((k-1)*NBy+l) + mfEy((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Ey_xz
          Ez_xz2(:,zz_cnt)=Ez_xz2(:,zz_cnt)+( pfEz((k-1)*NBy+l) + mfEz((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Ex_xz    
          
          Hx_xz2(:,zz_cnt)=Hx_xz2(:,zz_cnt)+( pfHx((k-1)*NBy+l) + mfHx((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Hx_xz
          Hy_xz2(:,zz_cnt)=Hy_xz2(:,zz_cnt)+( pfHy((k-1)*NBy+l) + mfHy((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Hy_xz
          Hz_xz2(:,zz_cnt)=Hz_xz2(:,zz_cnt)+( pfHz((k-1)*NBy+l) + mfHz((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % Hx_xz

      end; % for l
  end; % for k
      
  end; % zz_cnt
  
 
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

