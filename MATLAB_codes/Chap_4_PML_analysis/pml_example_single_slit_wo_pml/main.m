% FMM standard
% 2009 12 25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FMM with the extended S-matrix method
% written by Hwi Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. semi-infinite homogeneous space        - grating -  semi-infinite homogeneous space
% 2. semi-infinite homogeneous space        - grating -  semi-infinite inhomogeneous waveguide
% 3. semi-infinite inhomogeneous waveguide  - grating -  semi-infinite homogeneous space
% 4. semi-infinite inhomogeneous waveguide  - grating -  semi-infinite inhomogeneous waveguide

%% part 1-1
clear all;
close all;
clc;
tic;

addpath([pwd '\PRCWA_COM']);
%% part 1-2

% length unit

global nm;  % nano
global um;  % micro
global mm;  % mili

global k0;                                  % wavenumber
global c0; global w0;
global eps0; global mu0;


% zero-thickness buffer refractive index, permittivity, permeability
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

nano=10^-9;
micro=10^-6;
um=micro;
nm=nano;
%%
lambda=532*nm;          % operation wavelength
ni=1.5;                 % refractive index in region I
nf=1;                   % refractive index in region II
Tx=20.000001*micro;		% x-size of the computation cell
Ty=Tx;                  % y-size of the computation cell
nx=100;                 % x direction truncation order
ny=0;                   % y direction truncation order
theta=0*pi/180;         % incident angle
phi=0;                  % azimuthal angle
psi=0;                  % polarizatoin angle : TM=0, TE=pi/2



%% part 2

PRCWA_basic;                                %   3D structure
PRCWA_Gen_K;                                %   zero-thickness buffer

%% part 3 : structure modeling
PRCWA_Gen_PML2D;                            %   PML build
PRCWA_Gen_SPP_beaming;                     %   SPP beaming
%PRCWA_Gen_Y_branch;                         %   SPP Y branch

%% part 4-1 Layer S-matrix computation of layers

L=NBx*NBy;

Ta=zeros(2*L,2*L,Nlay); % left to rignt
Ra=zeros(2*L,2*L,Nlay); % left to right
Tb=zeros(2*L,2*L,Nlay); % right to left
Rb=zeros(2*L,2*L,Nlay); % right to left
Ca=zeros(4*L,2*L,Nlay); % left to right
Cb=zeros(4*L,2*L,Nlay); % left to right
tCa=zeros(4*L,2*L,Nlay); % left to right
tCb=zeros(4*L,2*L,Nlay); % left to right


%Off_diagonal_tensor_SMM;
Diagonal_SMM;


%--------------------------------------------------------------------------


%% part 4-2 S-matrix method of grating body


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

    end % for k

    tCa(:,:,laynt)=Ca(:,:,laynt)*inv(I-R_temp1b*R_temp2a)*T_temp1a;
    tCb(:,:,laynt)=Cb(:,:,laynt)+Ca(:,:,laynt)*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b;


    T_temp1a=TTa;
    R_temp1a=RRa;
    T_temp1b=TTb;
    R_temp1b=RRb;

    Ca=tCa;
    Cb=tCb;
end % laynt

TTa=T_temp1a;  % left-to-right transmission operator
RRa=R_temp1a;  % left-to-right reflection operator
TTb=T_temp1b;  % right-to-left transmission operator
RRb=R_temp1b;  % right-to-left reflection operator

%% part 4-2 S-matrix method of boundary matching

if Nlay==1
    Bdr_Smat_case1;      % 1. homogeneous space          - grating - homogeneous space
else
    Bdr_Smat_case2;       % 2. homogeneous space          - grating - inhomogeneous waveguide
end
%Bdr_Smat_case3;      % 3. inhomogeneous waveguide    - grating - homogeneous space
%Bdr_Smat_case4;      % 4. inhomogeneous waveguide    - grating - inhomogeneous waveguide

%% part 4-3 Field distribution and physical interpretation

if Nlay==1
    Field_visualization_3D_xz_case1_Lfree_Rfree_leftright;
else
    Field_visualization_3D_xz_case2_Lfree_Rwg_leftright;
end

% Field_visualization_3D_xz_case2_Lfree_Rwg_rightleft;
%
% Field_visualization_3D_yz_case2_Lfree_Rwg_leftright;
% Field_visualization_3D_yz_case2_Lfree_Rwg_rightleft;
%%
% save('beaming_example.mat','nx','pfEx','Tx','lambda');
% load('beaming_example.mat');
figure(2);plot((-nx:nx),abs(pfEx)/max(abs(pfEx)),'k');hold on;
set(gca,'fontsize',16);set(gca,'fontname','times new roman');
xlabel('Diffraction order');ylabel('|T_x| (arb. unit)');
axis([-nx nx 0 2]);
mm=floor(Tx/lambda);
line( mm*ones(1,10),linspace(0,2,10),'linestyle','--','color','k');
line(-mm*ones(1,10),linspace(0,2,10),'linestyle','--','color','k');

%%
toc


