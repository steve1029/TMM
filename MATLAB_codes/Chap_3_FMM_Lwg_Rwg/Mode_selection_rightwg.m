% Mode selection of the left waveguide

xx=3*[-Tx/2:Tx*0.002:Tx/2]';
z_start=0*um;
z_end=10*um;
z_inc=50*nano;
zz=(z_start:z_inc:z_end);
zn=length(zz);
y=0;
% region 1
laynt=Nlay; % left waveguide
    
pevalue=Peigvalue(1,:,laynt);
mevalue=Meigvalue(1,:,laynt);


[p_ind1 p_ind2]=sort( -real(pevalue) );
[m_ind1 m_ind2]=sort(  real(mevalue) );



pcnt=length(pevalue); % number of positive modes
mcnt=length(mevalue); % number of negative modes

zm=z_start;
zp=z_end;
       
zzp=zz-zm;
zzm=zz-zp;
    
Cp=zeros(pcnt,1);  % positive eigenmode
Cm=zeros(mcnt,1);  % negative eigenmode


modekind = 'positive';
%modekind = 'negative';
 
 
 switch lower(modekind)

  case 'positive'
% positive eigenmode  cp_cnt=1 (17)->1st symmetric TE, 
%cp_cnt=4 (14) -> 1st antisymmetric TE
for cp_cnt=1:20
    
    Cp=zeros(pcnt,1);
    Cp(p_ind2(cp_cnt))=1;
    
  pEy_xz=zeros(length(xx),length(zz));
  pEx_xz=zeros(length(xx),length(zz));
  pEz_xz=zeros(length(xx),length(zz));
  pHy_xz=zeros(length(xx),length(zz));
  pHx_xz=zeros(length(xx),length(zz));
  pHz_xz=zeros(length(xx),length(zz));     
  

    
  for zz_cnt=1:length(zzp)
   
    for k=1:2*L
        X(k,k)   =exp(pevalue(k)*zzp(zz_cnt));
    end;
    
    tmp1=X*Cp;
    pfEx=Pf_Ex(:,:,laynt)*tmp1;
    pfEy=Pf_Ey(:,:,laynt)*tmp1;
    pfEz=Pf_Ez(:,:,laynt)*tmp1;
    pfHx=Pf_Hx(:,:,laynt)*tmp1;
    pfHy=Pf_Hy(:,:,laynt)*tmp1;
    pfHz=Pf_Hz(:,:,laynt)*tmp1;
    
   

     for k=1:NBx
      for l=1:NBy
          
          pEx_xz(:,zz_cnt)=pEx_xz(:,zz_cnt)+( pfEx((k-1)*NBy+l))*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % pEx_xz    
          pEy_xz(:,zz_cnt)=pEy_xz(:,zz_cnt)+( pfEy((k-1)*NBy+l))*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % pEy_xz
          pEz_xz(:,zz_cnt)=pEz_xz(:,zz_cnt)+( pfEz((k-1)*NBy+l))*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % pEx_xz    
          
          pHx_xz(:,zz_cnt)=pHx_xz(:,zz_cnt)+( pfHx((k-1)*NBy+l))*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % pHx_xz
          pHy_xz(:,zz_cnt)=pHy_xz(:,zz_cnt)+( pfHy((k-1)*NBy+l))*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % pHy_xz
          pHz_xz(:,zz_cnt)=pHz_xz(:,zz_cnt)+( pfHz((k-1)*NBy+l))*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % pHx_xz

      end; % for l
  end; % for k
      
  end; % zz_cnt
  
  
figure(1); subplot(2,3,1);  imagesc(zz/um,xx/um,real(pEx_xz)); set(gca,'fontsize',16); colorbar; title('pEx-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);
           subplot(2,3,2);  imagesc(zz/um,xx/um,real(pEy_xz)); set(gca,'fontsize',16); colorbar; title('pEy-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);
           subplot(2,3,3);  imagesc(zz/um,xx/um,real(pEz_xz)); set(gca,'fontsize',16); colorbar; title('pEz-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);
           subplot(2,3,4);  imagesc(zz/um,xx/um,real(pHx_xz)); set(gca,'fontsize',16); colorbar; title('pHx-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);
           subplot(2,3,5);  imagesc(zz/um,xx/um,real(pHy_xz)); set(gca,'fontsize',16); colorbar; title('pHy-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);
           subplot(2,3,6);  imagesc(zz/um,xx/um,real(pHz_xz)); set(gca,'fontsize',16); colorbar; title('pHz-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);

  pause;
end; % for cp_cnt



 case 'negative'
     
% negative eigenmode
for cm_cnt=1:20

    Cm=zeros(mcnt,1);
    Cm(m_ind2(cm_cnt))=1;
    
  mEy_xz=zeros(length(xx),length(zz));
  mEx_xz=zeros(length(xx),length(zz));
  mEz_xz=zeros(length(xx),length(zz));
  mHy_xz=zeros(length(xx),length(zz));
  mHx_xz=zeros(length(xx),length(zz));
  mHz_xz=zeros(length(xx),length(zz));     

  for zz_cnt=1:length(zzp)
   
    for k=1:2*L
      Xinv(k,k)=exp(mevalue(k)*zzm(zz_cnt)); % there is no reflection mode
    end;
    
       
    tmp2=Xinv*Cm;
    mfEx=Mf_Ex(:,:,laynt)*tmp2;
    mfEy=Mf_Ey(:,:,laynt)*tmp2;
    mfEz=Mf_Ez(:,:,laynt)*tmp2;
    mfHx=Mf_Hx(:,:,laynt)*tmp2;
    mfHy=Mf_Hy(:,:,laynt)*tmp2;
    mfHz=Mf_Hz(:,:,laynt)*tmp2;
    
    
     for k=1:NBx
      for l=1:NBy
          
          mEx_xz(:,zz_cnt)=mEx_xz(:,zz_cnt)+( mfEx((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % mEx_xz    
          mEy_xz(:,zz_cnt)=mEy_xz(:,zz_cnt)+( mfEy((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % mEy_xz
          mEz_xz(:,zz_cnt)=mEz_xz(:,zz_cnt)+( mfEz((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % mEx_xz    
          
          mHx_xz(:,zz_cnt)=mHx_xz(:,zz_cnt)+( mfHx((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % mHx_xz
          mHy_xz(:,zz_cnt)=mHy_xz(:,zz_cnt)+( mfHy((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % mHy_xz
          mHz_xz(:,zz_cnt)=mHz_xz(:,zz_cnt)+( mfHz((k-1)*NBy+l) )*exp(j*( kx_vc(k)*xx+ky_vc(l)*y));   % mHx_xz

          end; % for l
      end; % for k
      
  end; % zz_cnt
 
  
figure(2); subplot(2,3,1);  imagesc(zz/um,xx/um,real(mEx_xz)); set(gca,'fontsize',16); colorbar; title('mEx-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);
           subplot(2,3,2);  imagesc(zz/um,xx/um,real(mEy_xz)); set(gca,'fontsize',16); colorbar; title('mEy-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);
           subplot(2,3,3);  imagesc(zz/um,xx/um,real(mEz_xz)); set(gca,'fontsize',16); colorbar; title('mEz-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);
           subplot(2,3,4);  imagesc(zz/um,xx/um,real(mHx_xz)); set(gca,'fontsize',16); colorbar; title('mHx-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);
           subplot(2,3,5);  imagesc(zz/um,xx/um,real(mHy_xz)); set(gca,'fontsize',16); colorbar; title('mHy-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);
           subplot(2,3,6);  imagesc(zz/um,xx/um,real(mHz_xz)); set(gca,'fontsize',16); colorbar; title('mHz-xz'); xlabel('z-axis(um)');ylabel('x-axis(um)'); colormap(hot); set(gca,'fontsize',16);

  pause;
end; % for cm_cnt


 end; % switch


