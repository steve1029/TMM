% Bdr_SMat_infr_outfr
function [Lfree_Tf,Lfree_Rb,Lfree_Tb,Lfree_Rf,Rfree_Tf,Rfree_Rb,Rfree_Tb,Rfree_Rf]=Bdr_SMat_infr_outfr()


global k0;                                  % wavenumber
global c0; global w0;
global eps0; global mu0;

global nano; global micro; global lambda; 
global n0; global epr0; global mur0; 
global ni; global epri; global muri;
global nf; global eprf; global murf;

global Tx; global Ty; global nx; global ny; 
global NBx; global NBy; global num_hx; global num_hy; global k0;
global kx_vc; global ky_vc; global kz_vc;

% input output free space
global kix; global kiy; global kiz; global kfz;
global kx_ref; global ky_ref; global kz_ref;
global kx_tra; global ky_tra; global kz_tra;



L=NBx*NBy;
% Buffere ???? ????
WW=zeros(2*L,2*L,2);	 	%W_L ... W_1
VV=zeros(2*L,2*L,2);    		%V_L ... V_1


Iden=zeros(2*L,2*L);
for k=1:2*L
Iden(k,k)=1;   
end;

% buffer
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


Wh=Iden;
Vh=KII/(w0*mu0);

% input freespace
Kin=zeros(2*L,2*L); % Kin matrix
for k=1:NBx
   for l=1:NBy
od_ind1=(k-1)*NBy+l;      
K11(od_ind1,od_ind1)=(kx_ref(k)*ky_ref(l))/kz_ref(k,l);
K12(od_ind1,od_ind1)=(kz_ref(k,l)^2+kx_ref(k)^2)/kz_ref(k,l);
K21(od_ind1,od_ind1)=-(ky_ref(l)^2+kz_ref(k,l)^2)/kz_ref(k,l);
K22(od_ind1,od_ind1)=-( ky_ref(l)*kx_ref(k) )/kz_ref(k,l);

   end;
end;

Kin(1:L, 1:L)=K11;
Kin(1:L, L+1:2*L)=K12;
Kin(L+1:2*L, 1:L)=K21;
Kin(L+1:2*L,L+1:2*L )=K22;

Win=Iden;
Vin=Kin/(w0*mu0);

% output freespace
Kout=zeros(2*L,2*L); % KII matrix
for k=1:NBx
   for l=1:NBy
od_ind1=(k-1)*NBy+l;      
K11(od_ind1,od_ind1)=(kx_tra(k)*ky_tra(l))/kz_tra(k,l);
K12(od_ind1,od_ind1)=(kz_tra(k,l)^2+kx_tra(k)^2)/kz_tra(k,l);
K21(od_ind1,od_ind1)=-(ky_tra(l)^2+kz_tra(k,l)^2)/kz_tra(k,l);
K22(od_ind1,od_ind1)=-( ky_tra(l)*kx_tra(k) )/kz_tra(k,l);

   end;
end;

Kout(1:L, 1:L)=K11;
Kout(1:L, L+1:2*L)=K12;
Kout(L+1:2*L, 1:L)=K21;
Kout(L+1:2*L,L+1:2*L )=K22;

clear K11;
clear K12;
clear K21;
clear K22;

Wout=Iden;
Vout=Kout/(w0*mu0);


WW(:,:,1)=Iden;
VV(:,:,1)=Vin;

WW(:,:,2)=Iden;
VV(:,:,2)=Vout;


% Input S-matrix
W0=WW(:,:,1);
V0=VV(:,:,1);
aa=inv(inv(W0)*Wh+inv(V0)*Vh);
bb=inv(inv(Wh)*W0+inv(Vh)*V0);

Lfree_Tf=2*aa;
Lfree_Rb=bb*(-inv(Wh)*W0+inv(Vh)*V0);
Lfree_Tb=2*bb;
Lfree_Rf=aa*(-inv(W0)*Wh+inv(V0)*Vh);

% Output S-matrix
Wf=WW(:,:,2);
Vf=VV(:,:,2);
aa=inv(inv(Wh)*Wf+inv(Vh)*Vf);
bb=inv(inv(Wf)*Wh+inv(Vf)*Vh);

Rfree_Tf=2*aa;
Rfree_Rb=bb*(-inv(Wf)*Wh+inv(Vf)*Vh);
Rfree_Tb=2*bb;
Rfree_Rf=aa*(-inv(Wh)*Wf+inv(Vh)*Vf);


%y=[inTf	inRb	inTb	inRf	outTf outRb	outTb	outRf];









