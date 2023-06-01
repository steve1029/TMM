%PRCWA_basic.m

% global vaiables-----------------------------------------------------
nano=10^-9;
micro=10^-6;
um=micro;
nm=nano;
    
% zero-thickness buffer
n0=1;                       % refractive index in zero-thickness buffer
epr0=n0^2;					% permittivity in zero-thickness buffer
mur0=1;                     % permeability in zero-thickness buffer

% freespace - grating - freespace
ni=1;                       %  refractive index in region I
epri=ni^2;                  % permittivity in region I
muri=1;                     % permeability in region I

nf=1;       				% refractive index in region II
eprf=nf^2;					% permittivity in region II
murf=1;						% permeability in region II

% photonic crystal structure
a=0.58*um;
Tx=15*a;
Ty=Tx;
Tz=a;

% basic constant
lambda=a/0.41;       % operating wavelength
c0=2.99792458*10^8;  % light speed
w0=2*pi*c0/lambda;   % angular frequency
eps0=1/(36*pi)*1e-9; % permittivity in free space
mu0=4*pi*10^-7;  

% for 2D structure
nx=14; 	 			% x direction truncation order
ny=0;     			% y direction truncation order
nz=14;               % z direction truncation order

% for 3D structure
% nx=1;              % x direction truncation order
% ny=1;              % y direction truncation order
% nz=1;              % z direction truncation order

NBx=2*nx+1; 		% number of supported orders of x axis
NBy=2*ny+1; 		% number of supported orders of y axis
NBz=2*nz+1;         % number of supported orders of z axis

num_hx=2*NBx-1; 
num_hy=2*NBy-1;           
num_hz=2*NBz-1;        

L=NBx*NBy;

%%%%%%%%%%%%%% Very important
aTx=Tx;
aTy=Ty;
aTz=Tz;

bTx=Tz;
bTy=Ty;
bTz=Tx;


