% PRCWA_Gen_K

theta=0*pi/180;         % incident angle
phi=0;      		    % azimuthal angle
psi=pi/2;               % polarization angle : 

Ux=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);  % incident wave�� Ex
Uy=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);  % incident wave�� Ey
Uz=-cos(psi)*sin(theta); 


%%% Incident E-field wavenumber %%%
k0=2*pi/lambda;                 % reference wavenumber
kx0=k0*ni*sin(theta)*cos(phi); % incident wavevector kx
ky0=k0*ni*sin(theta)*sin(phi); % incident wavevector ky
kz0=(k0^2-kx0^2-ky0^2)^0.5;          % incident wavevector kz

% harmonic max order num_hx, num_hy
% region I : incident & reflection 
kx_vc=zeros(NBx,1);     % supported kx modes in x direction: kxm
ky_vc=zeros(NBy,1);     % supported ky modes in y direction: kxn
kz_vc=zeros(NBx,NBy);   % supported kz modes in z direction: kzmn
R=zeros(NBx*NBy,1);     

for k=1:NBx
   for l=1:NBy  
	kx_vc(k)=kx0+(k-(nx+1))*(2*pi/Tx); 
	ky_vc(l)=ky0+(l-(ny+1))*(2*pi/Ty); 
	kz_vc(k,l)=k0*(n0^2-(kx_vc(k)/k0)^2-(ky_vc(l)/k0)^2)^0.5;
 end;
end;

% For A mode
Akx_vc=zeros(NBx,1);  % supported mode number in x direction: kxm
Aky_vc=zeros(NBy,1);  % supported mode number in y direction: kxn
Akz_vc=zeros(NBx,NBy);  % supported mode number in z direction: kzmn

for k=1:NBx
for l=1:NBy  
	Akx_vc(k)=(k-(nx+1))*(2*pi/aTx); 
	Aky_vc(l)=(l-(ny+1))*(2*pi/aTy); 
	Akz_vc(k,l)=k0*(n0^2-(Akx_vc(k)/k0)^2-(Aky_vc(l)/k0)^2)^0.5;	
end;
end;


% For B mode
Bkz_vc=zeros(NBx,1);  % supported mode number in x direction: kxm
Bky_vc=zeros(NBy,1);  % supported mode number in y direction: kxn
Bkx_vc=zeros(NBx,NBy);  % supported mode number in z direction: kzmn

for k=1:NBx
for l=1:NBy  
	Bkz_vc(k)=(k-(nx+1))*(2*pi/bTz); 
	Bky_vc(l)=(l-(ny+1))*(2*pi/bTy); 
	Bkx_vc(k,l)=k0*(n0^2-(Bkz_vc(k)/k0)^2-(Bky_vc(l)/k0)^2)^0.5;
end;
end;


Iden=zeros(2*L,2*L);
for k=1:2*L
Iden(k,k)=1;   
end;
      
% L+1 layer : region II
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




