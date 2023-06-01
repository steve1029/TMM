% positive eigenmode field visualization 
function [Ex_xz Ey_xz Ez_xz Hx_xz Hy_xz Hz_xz]=pEig_xz(pfEx,pfEy,pfEz,pfHx,pfHy,pfHz,pevalue,xx,zz)

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
          
          Ex_xz=Ex_xz+pfEx((k-1)*NBy+l)   *exp(j*( kx_vc(k)*x+ky_vc(l)*y)).*exp(pevalue*z);   % Ex_xz    
          Ey_xz=Ey_xz+pfEy((k-1)*NBy+l)   *exp(j*( kx_vc(k)*x+ky_vc(l)*y)).*exp(pevalue*z);   % Ey_xz
          Ez_xz=Ez_xz+pfEz((k-1)*NBy+l)   *exp(j*( kx_vc(k)*x+ky_vc(l)*y)).*exp(pevalue*z);   % Ex_xz    
          
          Hx_xz=Hx_xz+pfHx((k-1)*NBy+l)   *exp(j*( kx_vc(k)*x+ky_vc(l)*y)).*exp(pevalue*z);   % Hx_xz
          Hy_xz=Hy_xz+pfHy((k-1)*NBy+l)   *exp(j*( kx_vc(k)*x+ky_vc(l)*y)).*exp(pevalue*z);   % Hy_xz
          Hz_xz=Hz_xz+pfHz((k-1)*NBy+l)   *exp(j*( kx_vc(k)*x+ky_vc(l)*y)).*exp(pevalue*z);   % Hx_xz

      end; % for l
  end; % for k