% FMM single block structures written by H. Kim




%% STEP 1 
% All variables in workspace are cleard. All figures are closed and command
% window is initialized.

clear all;
close all;
clc;

% The library directories PRCWA_COM, FIELD_VISUAL, and STRUCTURE are added
% to path. The main directory of the code is also added to path.

addpath([pwd '\PRCWA_COM']);
addpath([pwd '\FIELD_VISUAL']);
addpath([pwd '\STRUCTURE']);

% Several basic variables are declared to global variabls. All global
% variables are set in PRCWA_basic.m, PRCWA_Gen_K.m in the subdirectory
% PRCWA_COM

% length unit
global nm;  % nano
global um;  % micro
global mm;  % mili

global k0;                                  % wavenumber
global c0; global w0;
global eps0; global mu0;
%ero-thickness buffer refractive index, permittivity, permeability 
global n0; global epr0; global mur0;
% refractive index, permittivity, permeability in free space I
global ni; global epri; global muri;
% refractive index, permittivity, permeability in free space II
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
lambda=532*nm; % Operating wavelength is set to 532nm

 % direct_=1 left-to-right characterization, 
 % direct_=2 right-to-left characterization
 direct_ = 1; % 1 = left-to-right  2 = right-to-left

 % The Fourier coefficients of the permittivity and peameability profiles of
% the target structure are calculated to Eqs. (4.1.1a) and (4.1.1b),
% respectively.


% In PRCWA_basic.m, the basic setting of Fourier harmonics orders, 
%supercell period, buffer layer material are done.
% Tx and Ty are the x-directional and y-directional supercell periods. The
% variables nx and ny are corresponding to the numbers M and N in Eqs.
% (4.1.1a) and (4.1.1b), respectively. In the example of this case, since
% nx is set to 0, the structure is the grating with spatial permittivity
% or permeability variations alont the y-axis.
% k0 is the wavenumber k0 is the wavenumber and kx0, ky0, and kz0 are the
% reference wavevector determined by the incidence angle theta and 
%azimuthal angle phi are the incidence angle and azimuthal angle of the 
%reference wavevector, respectively.

PRCWA_basic;                                %   3D structure

% The set of wavevectors in Eqs. (4.1.2a), (4.1.2b), and (4.1.2c) has to be
% prepared for the FMM. In PRCWA_Gen_K.m, kx_vc, ky_vc, and kz_vc are
% corresponding to Eqs. (4.1.2a), (4.1.2b) and (4.1.2c), respectively.

PRCWA_Gen_K;                                %   zero-thickness buffer

% In PRCWA_Gen_diagonal_BinaryGrating.m, the Fourier series coefficients of
% the structure are calculated. In Grating_gen_BinaryGrating_along_x.m, 
%rect_2D_mesh.m is the main function of the analytic Fourier series 
%coefficient of rectangle function. 
%Open PRCWA_Gen_BinaryGrating.m and see the annotation


%In the latter part of PRCWA_Gen_diagonal_BinaryGrating.m, 
%the calculated Fourier series are saved as the form of tensor

%According to the tensor structure, two forms are distinguished 
%in the FMM package.In PRCWA_Gen_offdiagonal_BinaryGrating.m, the full 
%off-diagonal tensor terms are to have appropriate values. 
%In the code example, the isotropic material grating is exemplified and so 
%all off-diagonal terms are set to 0.

PRCWA_Gen_BinaryGrating;    %  full off-diagonal anisotropic material       

%% STEP 2 Block S-matrix computation of single block structures

L=NBx*NBy;

TTa=zeros(2*L,2*L,Nlay); % left to rignt transmission S-matrix
RRa=zeros(2*L,2*L,Nlay); % left to right reflection S-matrix
TTb=zeros(2*L,2*L,Nlay); % right to left transmission S-matrix
RRb=zeros(2*L,2*L,Nlay); % right to left reflection S-matrix
Ca=zeros(4*L,2*L,Nlay);  % left to right coupling coefficient matrix operator
Cb=zeros(4*L,2*L,Nlay);  % left to right coupling coefficient matrix operator

laynt=1;

% With the structural modeling and discete wavevector set, we can analyze 
% the internal electromagnetic Bloch eigenmodes by solving numerical 
% eigenvalue equation. As described in Chapter 3, the Fourier modal 
% representation of the eigenmodes are described as the eigensolution 
% of the Maxwell eigenvalue equations.

% Two MATLAB function for the analysis of single block structures
% (i) FMM_single_block.m
% (ii) FMM_single_block_tensor.m
% Both functions generate the S-matrix, coupling coefficient operator, 
% pseudo-Fourier series coefficients of eigenmodes, and eigenvalues with 
% the input of the permittivity and permeability profiles.
% In the example, FMM_single_block.m is used with PRCWA_Gen_diagonal_BinaryGrating.m 
% and FMM_single_block_tensor.m is used with
% PRCWA_Gen_offdiagonal_BinaryGrating.m.

% Open two m-files and see the annotation

Diagonal_SMM;
%Off_diagonal_tensor_SMM;


%% STEP3 Field visualization

L=NBx*NBy;

K=zeros(2*L,2*L); % KII matrix
for k=1:NBx
   for l=1:NBy
od_ind1=(k-1)*NBy+l;      
K11(od_ind1,od_ind1)=(kx_vc(k)*ky_vc(l))/kz_vc(k,l);
K12(od_ind1,od_ind1)=(kz_vc(k,l)^2+kx_vc(k)^2)/kz_vc(k,l);
K21(od_ind1,od_ind1)=-(ky_vc(l)^2+kz_vc(k,l)^2)/kz_vc(k,l);
K22(od_ind1,od_ind1)=-( ky_vc(l)*kx_vc(k) )/kz_vc(k,l);

   end;
end;

K(1:L, 1:L)=K11;
K(1:L, L+1:2*L)=K12;
K(L+1:2*L, 1:L)=K21;
K(L+1:2*L,L+1:2*L )=K22;

Iden=zeros(2*L,2*L);
for k=1:2*L
Iden(k,k)=1;   
end;

Wh=Iden;
Vh=K/(w0*mu0);

% The field visualization of six E & M components is basic task.
% We prepared five field visualization MATLAB codes to test the S-matrix
% formulation and provide the reference MATLAB codes.

switch direct_

    case 1      % left-to-right
        
        % polarizatoin angle : TM=0, TE=pi/2
psi=0;      			
tm_Ux=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);  
tm_Uy=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);  
tm_Uz=-cos(psi)*sin(theta);  

psi=pi/2;
te_Ux=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);  
te_Uy=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);  
te_Uz=-cos(psi)*sin(theta); 
        
Field_visualization_3D_xz_leftright;
%Field_visualization_3D_yz_leftright;
        
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
            
Field_visualization_3D_xz_rightleft;
%Field_visualization_3D_yz_rightleft;
         
        
end;

%Field_visualization_3D_xy;

