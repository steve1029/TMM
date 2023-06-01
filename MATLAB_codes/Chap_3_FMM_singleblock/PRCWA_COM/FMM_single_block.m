% FMM_single_block for off-diagonal anisotropic block
% Bidirectional characterization

% Eqs. (4.1.4a) and (4.1.5b) are corresponding to SA and SB, respectively.
% Eq. (4.1.5c) is implemented by St in MATLAB code. The Maxwell eigenvalue
% equation (4.1.5c) is numerically solved and the resulting eigenvalues and
% eigenvectors are saved to W and Dt, respectively. This part is described
% in MATLALB code lines 103~159. Lalanne's experience rule of use of
% reciprocal permittivity and permeablity is employed.
% On the other hand, when the material is fully anisotropic,
% use FMM_single_block_tensor.m. Open the code and see the annotation.


function [Ta,Ra,Tb,Rb,CCa,CCb,pfEx,pfEy,pfEz,pfHx,pfHy,pfHz,pevalue,mfEx,mfEy,mfEz,mfHx,mfHy,mfHz,mevalue]=FMM_single_block(lay_thick,str_tensor,alpha_tm,beta_tm)

global c0; global w0;
global eps0; global mu0;

global nano; global micro; global lambda; 
global n0; global epr0; global mur0;
global eps0; global mu0;
global Tx; global Ty; global nx; global ny; 
global NBx; global NBy; global num_hx; global num_hy; global k0;
global kx_vc; global ky_vc; global kz_vc;


L=NBx*NBy; 

Teigvalue=zeros(1,2*L);

Kx=zeros(L,L);
Ky=zeros(L,L);
for k=1:NBx 
   for l=1:NBy  
      od_ind=(k-1)*NBy+l;
      Kx(od_ind,od_ind)=kx_vc(k)/k0;
      Ky(od_ind,od_ind)=ky_vc(l)/k0;   
   end;
end;


% Toeplitz matrix
Exx=zeros(L,L);	 
Eyy=zeros(L,L);
Ezz=zeros(L,L);

Axx=zeros(L,L);   
Ayy=zeros(L,L);
Azz=zeros(L,L);

Gxx=zeros(L,L);
Gyy=zeros(L,L);
Gzz=zeros(L,L);

Bxx=zeros(L,L);
Byy=zeros(L,L);
Bzz=zeros(L,L);

eps_xx=str_tensor.eps_xx;
eps_yy=str_tensor.eps_yy;
eps_zz=str_tensor.eps_zz;

aps_xx=str_tensor.aps_xx;
aps_yy=str_tensor.aps_yy;
aps_zz=str_tensor.aps_zz;

mu_xx=str_tensor.mu_xx;
mu_yy=str_tensor.mu_yy;
mu_zz=str_tensor.mu_zz;

bu_xx=str_tensor.bu_xx;
bu_yy=str_tensor.bu_yy;
bu_zz=str_tensor.bu_zz;


for k=1:NBx
   for l=1:NBy
      od_ind1=(k-1)*NBy+l;
      for kk=1:NBx
      	for ll=1:NBy   
            od_ind2=(kk-1)*NBy+ll;		% od_ind2-od_ind1=
            
            % permittivity
            Exx(od_ind1,od_ind2)=eps_xx(k-kk+NBx,l-ll+NBy); 	
            Eyy(od_ind1,od_ind2)=eps_yy(k-kk+NBx,l-ll+NBy); 
            Ezz(od_ind1,od_ind2)=eps_zz(k-kk+NBx,l-ll+NBy); 
            
         
            Axx(od_ind1,od_ind2)=aps_xx(k-kk+NBx,l-ll+NBy); 	
            Ayy(od_ind1,od_ind2)=aps_yy(k-kk+NBx,l-ll+NBy); 
            Azz(od_ind1,od_ind2)=aps_zz(k-kk+NBx,l-ll+NBy); 

                       
            % permeability
			Gxx(od_ind1,od_ind2)=mu_xx(k-kk+NBx,l-ll+NBy); 
            Gyy(od_ind1,od_ind2)=mu_yy(k-kk+NBx,l-ll+NBy); 
		    Gzz(od_ind1,od_ind2)=mu_zz(k-kk+NBx,l-ll+NBy);     
            
            
            Bxx(od_ind1,od_ind2)=bu_xx(k-kk+NBx,l-ll+NBy); 
            Byy(od_ind1,od_ind2)=bu_yy(k-kk+NBx,l-ll+NBy); 
	    	Bzz(od_ind1,od_ind2)=bu_zz(k-kk+NBx,l-ll+NBy);         

                      
      	end;
      end;
   end;
end;


alpha=alpha_tm; % convergence factor
beta=beta_tm;

AE_x=zeros(L,L);
EA_x=zeros(L,L);

AE_y=zeros(L,L);
EA_y=zeros(L,L);

AE_z=zeros(L,L);
EA_z=zeros(L,L);

AE_x=alpha*inv(Axx)+(1-alpha)*Exx;
EA_x=alpha*Exx+(1-alpha)*inv(Axx);

AE_y=alpha*inv(Ayy)+(1-alpha)*Eyy;
EA_y=alpha*Eyy+(1-alpha)*inv(Ayy);

AE_z=alpha*inv(Azz)+(1-alpha)*Ezz;
EA_z=alpha*Ezz+(1-alpha)*inv(Azz);

BG_x=zeros(L,L);
GB_x=zeros(L,L);

BG_y=zeros(L,L);
GB_y=zeros(L,L);

BG_z=zeros(L,L);
GB_z=zeros(L,L);


BG_x=beta*inv(Bxx)+(1-beta)*Gxx;
GB_x=beta*Gxx+(1-beta)*inv(Bxx);


BG_y=beta*inv(Byy)+(1-beta)*Gyy;
GB_y=beta*Gyy+(1-beta)*inv(Byy);


BG_z=beta*inv(Bzz)+(1-beta)*Gzz;
GB_z=beta*Gzz+(1-beta)*inv(Bzz);


% System Matrix

SA=zeros(2*L, 2*L);
SB=zeros(2*L,2*L);         




SA=[ (Ky)*inv(Ezz)*(Kx) 		 BG_x-(Ky)*inv(Ezz)*(Ky)
     (Kx)*inv(Ezz)*(Kx)-GB_y   -(Kx)*inv(Ezz)*(Ky)];

SB=[ (Ky)*inv(Gzz)*(Kx)  		 AE_x-(Ky)*inv(Gzz)*(Ky)
     (Kx)*inv(Gzz)*(Kx)-EA_x   -(Kx)*inv(Gzz)*(Ky)];

St=k0^2*SA*SB;
clear SB;

[W,Dt]=eig(St);  
clear St;


eig_value=zeros(1,2*L);
for k=1:2*L
eig_value(k)=Dt(k,k);   
end;

c_eig=eig_value;
for k=1:2*L
   eig_value(k)=c_eig(k)^0.5;  
   if real(eig_value(k))>0
       eig_value(k)=-eig_value(k);
   end
end;

mcnt=0;
pcnt=0;

%------------------------------------------------------------------------
for k=1:2*L
     
  % positive mode
       
       pcnt=pcnt+1;
       pevalue(pcnt)=eig_value(k);
       pW(:,pcnt)=W(:,k);
  % negative mode
   
       mcnt=mcnt+1;
       mevalue(mcnt)=-eig_value(k);
       mW(:,mcnt)=W(:,k);

end;

pQ=zeros(2*L,2*L);
mQ=zeros(2*L,2*L);

for k=1:2*L
  pQ(k,k)=pevalue(k);
  mQ(k,k)=mevalue(k);
end;


Z=inv(k0*SA);
pV=Z*pW*pQ;
mV=Z*mW*mQ;
clear Z;
clear SA;
 
% Fourier coefficients (pfEx,pfEy,pfEz,pfHx,pfHy,pfHz) , (mfEx,mfEy,mfEz,mfHx,mfHy,mfHz) 
 
    pfEy=pW(1:L,:);                                            %     pfEy
    pfEx=pW(L+1:2*L,:);                                        %     pfEx
    pfHy=pV(1:L,:);                                            %     pfHy
    pfHx=pV(L+1:2*L,:);                                        %     pfHx  
    pfEz=inv(Ezz)*(j*Ky*pfHx-j*Kx*pfHy);                       %     pfEz
    pfHz=inv(Gzz)*(j*Ky*pfEx-j*Kx*pfEy);                       %     pfHz

    pfHy=j*(eps0/mu0)^0.5*pfHy;                                %     pfHy
    pfHx=j*(eps0/mu0)^0.5*pfHx;                                %     pfHx
    pfHz=j*(eps0/mu0)^0.5*pfHz;                                %     pfHz
    

    mfEy=mW(1:L,:);                                            %     mfEy
    mfEx=mW(L+1:2*L,:);                                        %     mfEx
    mfHy=mV(1:L,:);                                            %     mfHy
    mfHx=mV(L+1:2*L,:);                                        %     mfHx
    mfEz=inv(Ezz)*(j*Ky*mfHx-j*Kx*mfHy);                       %     mfEz
    mfHz=inv(Gzz)*(j*Ky*mfEx-j*Kx*mfEy);                       %     mfHz
    
    mfHy=j*(eps0/mu0)^0.5*mfHy;                                %     mfHy
    mfHx=j*(eps0/mu0)^0.5*mfHx;                                %     mfHx
    mfHz=j*(eps0/mu0)^0.5*mfHz;    
    
%% single layer S-matrix
zm=0;
zp=lay_thick;

Wp_zp=zeros(2*L,pcnt); 
Wm_zp=zeros(2*L,mcnt);
Vp_zp=zeros(2*L,pcnt);
Vm_zp=zeros(2*L,mcnt);

Wp_zm=zeros(2*L,pcnt); 
Wm_zm=zeros(2*L,mcnt);
Vp_zm=zeros(2*L,pcnt);
Vm_zm=zeros(2*L,mcnt);


Wp_zm=sWp_gen(pW,pevalue,pcnt,L,zm-zm); Wm_zm=sWm_gen(mW,mevalue,mcnt,L,zm-zp);
Vp_zm=sVp_gen(pV,pevalue,pcnt,L,zm-zm); Vm_zm=sVm_gen(mV,mevalue,mcnt,L,zm-zp);

Wp_zp=sWp_gen(pW,pevalue,pcnt,L,zp-zm); Wm_zp=sWm_gen(mW,mevalue,mcnt,L,zp-zp);
Vp_zp=sVp_gen(pV,pevalue,pcnt,L,zp-zm); Vm_zp=sVm_gen(mV,mevalue,mcnt,L,zp-zp);



%--------------------------------------------------------------------------
Iden=zeros(2*L,2*L);
for k=1:2*L
Iden(k,k)=1;   
end;

% region II
KII=zeros(2*L,2*L); % KII matrix
for k=1:NBx
   for l=1:NBy
od_ind1=(k-1)*NBy+l;      
K11(od_ind1,od_ind1)=(kx_vc(k)*ky_vc(l))/kz_vc(k,l);
K12(od_ind1,od_ind1)=(kz_vc(k,l)^2+kx_vc(k)^2)/kz_vc(k,l);
K21(od_ind1,od_ind1)=-(ky_vc(l)^2+kz_vc(k,l)^2)/kz_vc(k,l);
K22(od_ind1,od_ind1)=-( ky_vc(l)*kx_vc(k) )/kz_vc(k,l);
   end;
end;

KII(1:L, 1:L)=K11;
KII(1:L, L+1:2*L)=K12;
KII(L+1:2*L, 1:L)=K21;
KII(L+1:2*L,L+1:2*L )=K22;
clear K11;
clear K12;
clear K21;
clear K22;


Wh=Iden;
Vh=KII/(w0*mu0);
%--------------------------------------------------------------------------



 % left-to-right
U=eye(2*L,2*L);

S11=inv(Wh)*Wp_zm+inv(Vh)*Vp_zm;
S12=inv(Wh)*Wm_zm+inv(Vh)*Vm_zm;
S21=inv(Wh)*Wp_zp-inv(Vh)*Vp_zp;
S22=inv(Wh)*Wm_zp-inv(Vh)*Vm_zp;

S=[S11 S12; S21 S22];


clear S11;
clear S12;
clear S21;
clear S22;
D=[2*U;zeros(2*L,2*L)];

CCa=inv(S)*D;
Cap=CCa(1:2*L,:);
Cam=CCa(2*L+1:4*L,:);
Ra=inv(Wh)*(Wp_zm*Cap+Wm_zm*Cam-Wh*U);
Ta=inv(Wh)*(Wp_zp*Cap+Wm_zp*Cam);

 % right-to-left

 S11=inv(Wh)*Wp_zm+inv(Vh)*Vp_zm;
 S12=inv(Wh)*Wm_zm+inv(Vh)*Vm_zm;
 S21=inv(Wh)*Wp_zp-inv(Vh)*Vp_zp;
 S22=inv(Wh)*Wm_zp-inv(Vh)*Vm_zp;

S=[S11 S12 ; S21 S22];
clear S11;
clear S12;
clear S21;
clear S22;

D=[zeros(2*L,2*L);2*U];
CCb=inv(S)*D;
Cbp=CCb(1:2*L,:);
Cbm=CCb(2*L+1:4*L,:);

Rb=inv(Wh)*(Wp_zp*Cbp+Wm_zp*Cbm-Wh*U);
Tb=inv(Wh)*(Wp_zm*Cbp+Wm_zm*Cbm);



