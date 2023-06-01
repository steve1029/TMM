% Gen_K_space

theta=0*pi/180;                 % incident angle
phi=0;                          % azimuthal angle

%%% Incident E-field wavenumber %%%
k0=2*pi/lambda;                 % reference wavenumber
kx0=k0*ni*sin(theta)*cos(phi);  % incident wavevector kx
ky0=k0*ni*sin(theta)*sin(phi);  % incident wavevector ky
kz0=(k0^2-kx0^2-ky0^2)^0.5;     % incident wavevector kz

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

L=NBx*NBy;

Iden=zeros(2*L,2*L);
for k=1:2*L
Iden(k,k)=1;   
end;



