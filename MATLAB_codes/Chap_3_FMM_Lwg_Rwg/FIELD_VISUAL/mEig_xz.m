% negative eigenmode field visualization
function [Ex_xz Ey_xz Ez_xz Hx_xz Hy_xz Hz_xz]=mEig_xz(mfEx,mfEy,mfEz,mfHx,mfHy,mfHz,mevalue,xx,zz)

global kx_vc; global ky_vc; global kz_vc;
global NBx; global NBy;                     
L=NBx*NBy;

  Ex_xz=zeros(length(xx),length(zz));
  Ey_xz=zeros(length(xx),length(zz));
  Ez_xz=zeros(length(xx),length(zz));
  
  Hx_xz=zeros(length(xx),length(zz));
  Hy_xz=zeros(length(xx),length(zz));
  Hz_xz=zeros(length(xx),length(zz));     
  
  [z x]=meshgrid(zz,xx); % grating region
  y=0;
  
  for k=1:NBx
      for l=1:NBy
          
          Ex_xz=Ex_xz+mfEx((k-1)*NBy+l)   *exp(j*( kx_vc(k)*x+ky_vc(l)*y )).*exp(mevalue*z);   % Ex_xz    
          Ey_xz=Ey_xz+mfEy((k-1)*NBy+l)   *exp(j*( kx_vc(k)*x+ky_vc(l)*y )).*exp(mevalue*z);   % Ey_xz
          Ez_xz=Ez_xz+mfEz((k-1)*NBy+l)   *exp(j*( kx_vc(k)*x+ky_vc(l)*y )).*exp(mevalue*z);   % Ex_xz    
          
          Hx_xz=Hx_xz+mfHx((k-1)*NBy+l)   *exp(j*( kx_vc(k)*x+ky_vc(l)*y )).*exp(mevalue*z);   % Hx_xz
          Hy_xz=Hy_xz+mfHy((k-1)*NBy+l)   *exp(j*( kx_vc(k)*x+ky_vc(l)*y )).*exp(mevalue*z);   % Hy_xz
          Hz_xz=Hz_xz+mfHz((k-1)*NBy+l)   *exp(j*( kx_vc(k)*x+ky_vc(l)*y )).*exp(mevalue*z);   % Hx_xz

      end; % for l
  end; % for k