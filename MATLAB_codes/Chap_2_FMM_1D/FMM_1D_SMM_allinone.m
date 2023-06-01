
% FMM standard code written by H. Kim
% semi-infinite homogeneous space I - multilayer -  semi-infinite homogeneous space II       

clear all;
close all;
clc;

% addpath([pwd '\PRCWA_COM']);
% addpath([pwd '\FIELD_VISUAL']);
% addpath([pwd '\STRUCTURE']);

%% STEP 1 : wavevector setting and structure modeling

% length unit
global nm;  % nano
global um;  % micro
global mm;  % mili

% basic parameters
global k0;               % wavenumber
global c0; global w0;    % speed of light, angular frequency
global eps0; global mu0; % vaccum permittivity & permeability

% refractive index, permittivity, permeability in zero-thickness buffer layer 
global n0; global epr0; global mur0;
% refractive index, permittivity, permeability in homogeneous space I
global ni; global epri; global muri;
% refractive index, permittivity, permeability in homogeneous space II
global nf; global eprf; global murf;        

% x-directional supercell period, y-directional supercell period
global Tx; global Ty;                       
global nx; global ny;                       
% # of x-direction Fourier harmonics, # of y-directional Fourier harmonics
global NBx; global NBy;                     
global num_hx; global num_hy;               
global kx_vc; global ky_vc; global kz_vc;

% input output free space
global kix; global kiy; global kiz; 
global kfz;
global kx_ref; global ky_ref; global kz_ref;
global kx_tra; global ky_tra; global kz_tra;

nm=1e-9;
lambda=532*nm;

 % direct_=1 left-to-right characterization, 
 % direct_=2 right-to-left characterization
 direct_ =1; % 1 = left-to-right , 2 = right-to-left

% 3D structure
% PRCWA_basic.m
% global vaiables-----------------------------------------------------
nano=10^-9;
micro=10^-6;
um=micro;
nm=nano;

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
ni=1.5;                       %  refractive index in region I
epri=ni^2;                  % permittivity in region I
muri=1;                     % permeability in region I

nf=1.7;       				% refractive index in region II
eprf=nf^2;					% permittivity in region II
murf=1;						% permeability in region II

Tx=5.000001*micro;			
Ty=Tx;


% zero-thickness buffer
% Gen_K_space

theta=10*pi/180;                % incident angle
phi=0;                          % azimuthal angle

%%% Incident E-field wavenumber %%%
switch direct_
    
    case 1      % left-to-right
        k0=2*pi/lambda;                     % reference wavenumber
        kx0=k0*ni*sin(theta)*cos(phi);      % incident wavevector kx
        ky0=k0*ni*sin(theta)*sin(phi);      % incident wavevector ky
        kz0=(k0^2-kx0^2-ky0^2)^0.5;         % incident wavevector kz

    case 2      % right-to-left
        k0=2*pi/lambda;                     % reference wavenumber
        kx0=k0*nf*sin(theta)*cos(phi);      % incident wavevector kx
        ky0=k0*nf*sin(theta)*sin(phi);      % incident wavevector ky
        kz0=(k0^2-kx0^2-ky0^2)^0.5;         % incident wavevector kz

end;
 
kx_vc=kx0; 
ky_vc=ky0; 
kz_vc=k0*(n0^2-(kx_vc/k0)^2-(ky_vc/k0)^2)^0.5;

L=1;
Iden=zeros(2*L,2*L);
for k=1:2*L
Iden(k,k)=1;   
end;


% The example structure is a multi-layer structure with randomly generated
% refractive indices
% SPP Y branch
% Layer structure ------------------------------------------
Nlay=1+0;

% material permittivity
grating_thick=5*um;            % grating region

lay_thick=zeros(Nlay,1);      %   layer 
lay_thick(1)=0;                 %   left semi-infinite waveguide
for laynt=2:Nlay-1
  lay_thick(laynt)=grating_thick/(Nlay-2);  
end;
lay_thick(Nlay)=0;                 %   right semi-infinite waveguide

ac_thick=zeros(1,Nlay); 
for cnt1=1:Nlay
	 cnt2=1;
    while cnt2 <= cnt1
       ac_thick(cnt1)=ac_thick(cnt1)+lay_thick(cnt2);
       cnt2=cnt2+1;
    end;
end; % for cnt1


% Permittivity & Permeability
eps_L=zeros(1,Nlay);
aps_L=zeros(1,Nlay);
mu_L =zeros(1,Nlay);
bu_L =zeros(1,Nlay);

eps_L(1,:)=(1+rand(1,Nlay)+0.0000001*i).^2; % random permittivity
aps_L=1./eps_L;
mu_L(1,:)=1;  % random permeability
bu_L=1./mu_L;

eps_L(1,1)   =(ni+0.000001*i)^2;
eps_L(1,Nlay)=(nf+0.000001*i)^2;

alpha_tm=0;
beta_tm=0;

%---------------------- e_h, a_h, g_h, b_h -----------------------
eps_xx=zeros(1,Nlay); eps_xy=zeros(1,Nlay); eps_xz=zeros(1,Nlay);   % permittivity
eps_yx=zeros(1,Nlay); eps_yy=zeros(1,Nlay); eps_yz=zeros(1,Nlay);
eps_zx=zeros(1,Nlay); eps_zy=zeros(1,Nlay); eps_zz=zeros(1,Nlay); 

aps_xx=zeros(1,Nlay); aps_xy=zeros(1,Nlay); aps_xz=zeros(1,Nlay);   % reciprocal permittivity 
aps_yx=zeros(1,Nlay); aps_yy=zeros(1,Nlay); aps_yz=zeros(1,Nlay);
aps_zx=zeros(1,Nlay); aps_zy=zeros(1,Nlay); aps_zz=zeros(1,Nlay); 

mu_xx=zeros(1,Nlay);  mu_xy=zeros(1,Nlay);  mu_xz=zeros(1,Nlay);    % permeability
mu_yx=zeros(1,Nlay);  mu_yy=zeros(1,Nlay);  mu_yz=zeros(1,Nlay);
mu_zx=zeros(1,Nlay);  mu_zy=zeros(1,Nlay);  mu_zz=zeros(1,Nlay); 

bu_xx=zeros(1,Nlay);  bu_xy=zeros(1,Nlay);  bu_xz=zeros(1,Nlay);    % reciprocal permeability
bu_yx=zeros(1,Nlay);  bu_yy=zeros(1,Nlay);  bu_yz=zeros(1,Nlay);
bu_zx=zeros(1,Nlay);  bu_zy=zeros(1,Nlay);  bu_zz=zeros(1,Nlay); 


for lnt=1:Nlay    %diagonal anisotropic material multilayer
    
       eps_xx(1,lnt)=eps_L(1,lnt);
       eps_yy(1,lnt)=eps_L(1,lnt);
       eps_zz(1,lnt)=eps_L(1,lnt);
       aps_xx(1,lnt)=aps_L(1,lnt);
       aps_yy(1,lnt)=aps_L(1,lnt);
       aps_zz(1,lnt)=aps_L(1,lnt);
       mu_xx(1,lnt)=  mu_L(1,lnt);
       mu_yy(1,lnt)=  mu_L(1,lnt);
       mu_zz(1,lnt)=  mu_L(1,lnt);
       bu_xx(1,lnt)=  bu_L(1,lnt);
       bu_yy(1,lnt)=  bu_L(1,lnt);
       bu_zz(1,lnt)=  bu_L(1,lnt);
               
end; % for lnt


%% STEP 2 Block S-matrix computation of single block structures

L=1;

Ta=zeros(2*L,2*L,Nlay); % left to rignt
Ra=zeros(2*L,2*L,Nlay); % left to right
Tb=zeros(2*L,2*L,Nlay); % right to left
Rb=zeros(2*L,2*L,Nlay); % right to left
Ca=zeros(4*L,2*L,Nlay); % left to right
Cb=zeros(4*L,2*L,Nlay); % right to left
tCa=zeros(4*L,2*L,Nlay); % left to right
tCb=zeros(4*L,2*L,Nlay); % right to left

%Diagonal_SMM;               % diagonal anisotropy
Off_diagonal_tensor_SMM;   % off-diagonal anisotropy  

%%  STEP3 S-matrix method 
% The obtained S-matrix and coupling coefficient matrix operator of single
% blocks are combined to generate the S-matrix and coupling coefficient
% operator of interconnected multi-block structures by the Redheffer star
% product...
     
I=eye(2*L,2*L);
T_temp1a=Ta(:,:,1);
R_temp1a=Ra(:,:,1);
T_temp1b=Tb(:,:,1);
R_temp1b=Rb(:,:,1);
    
 %% Important
   tCa=Ca;
   tCb=Cb;
   %%
   
   
for laynt=2:Nlay
   
   %laynt
   
   T_temp2a=Ta(:,:,laynt);
   R_temp2a=Ra(:,:,laynt);
   T_temp2b=Tb(:,:,laynt);
   R_temp2b=Rb(:,:,laynt);

	RRa=(R_temp1a+T_temp1b*inv(I-R_temp2a*R_temp1b)*R_temp2a*T_temp1a);
	TTa=T_temp2a*inv(I-R_temp1b*R_temp2a)*T_temp1a;

	RRb=(R_temp2b+T_temp2a*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b);
    TTb=T_temp1b*inv(I-R_temp2a*R_temp1b)*T_temp2b;
   
   
   for k=1:laynt-1
   	
   tCa(:,:,k)=Ca(:,:,k)+Cb(:,:,k)*inv(I-R_temp2a*R_temp1b)*R_temp2a*T_temp1a;
   tCb(:,:,k)=Cb(:,:,k)*inv(I-R_temp2a*R_temp1b)*T_temp2b;

   end; % for k
   
   tCa(:,:,laynt)=Ca(:,:,laynt)*inv(I-R_temp1b*R_temp2a)*T_temp1a;
   tCb(:,:,laynt)=Cb(:,:,laynt)+Ca(:,:,laynt)*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b;
   

   T_temp1a=TTa;
   R_temp1a=RRa;
   T_temp1b=TTb;
   R_temp1b=RRb;
   
   Ca=tCa;
   Cb=tCb; 
end; % laynt
     
 TTa=T_temp1a;  % left-to-right transmission operator
 RRa=R_temp1a;  % left-to-right reflection operator
 TTb=T_temp1b;  % right-to-left transmission operator
 RRb=R_temp1b;  % right-to-left reflection operator
 
%% STEP4 Half-infinite block interconnection & Field visualization
 
 switch direct_
    
    case 1      % left-to-right
            % polarizatoin angle : TM=0, TE=pi/2
            psi=0;      			
            tm_Ux=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);  % incident wave�� Ex
            tm_Uy=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);  % incident wave�� Ey
            tm_Uz=-cos(psi)*sin(theta);  

            psi=pi/2;
            te_Ux=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);  % incident wave�� Ex
            te_Uy=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);  % incident wave�� Ey
            te_Uz=-cos(psi)*sin(theta); 

            % Ux=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);  % incident wave�� Ex
            % Uy=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);  % incident wave�� Ey
            % Uz=-cos(psi)*sin(theta); 

    case 2      % right-to-left
                    % polarizatoin angle : TM=0, TE=pi/2
            psi=0;      			
            tm_Ux=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);  % incident wave�� Ex
            tm_Uy=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);  % incident wave�� Ey
            tm_Uz=cos(psi)*sin(theta);  

            psi=pi/2;
            te_Ux=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);  % incident wave�� Ex
            te_Uy=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);  % incident wave�� Ey
            te_Uz=cos(psi)*sin(theta); 

            % Ux=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);  % incident wave�� Ex
            % Uy=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);  % incident wave�� Ey
            % Uz=cos(psi)*sin(theta); 

end;

Bdr_Smat_case1;   % 1. homogeneous space          - grating - homogeneous space    