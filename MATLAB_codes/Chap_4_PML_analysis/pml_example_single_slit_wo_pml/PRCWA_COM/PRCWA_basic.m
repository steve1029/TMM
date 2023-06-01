%PRCWA_basic.m

% global vaiables-----------------------------------------------------
% nano=10^-9;
% micro=10^-6;
% um=micro;
% nm=nano;

% basic constant
c0=2.99792458*10^8; % light speed
w0=2*pi*c0/lambda;   % angular frequency
eps0=1/(36*pi)*1e-9; % permittivity in free space
mu0=4*pi*10^-7;      % permeability in free space

% zero-thickness buffer
n0=1;                       % refractive index in zero-thickness buffer
epr0=n0^2;					% permittivity in zero-thickness buffer
mur0=1;                     % permeability in zero-thickness buffer

% freespace - grating - freespace
% ni=1;                       %  refractive index in region I
epri=ni^2;                  % permittivity in region I
muri=1;                     % permeability in region I

% nf=1;       				% refractive index in region II
eprf=nf^2;					% permittivity in region II
murf=1;						% permeability in region II

% Tx=20.000001*micro;			
% Ty=Tx;

% for 3D structure
% nx=100;              % x direction truncation order
% ny=0;               % y direction truncation order
NBx=2*nx+1; 		% number of supported orders of x axis
NBy=2*ny+1; 		% number of supported orders of y axis
num_hx=2*NBx-1;     % M+1 <-> odd number
num_hy=2*NBy-1;     % odd number          
