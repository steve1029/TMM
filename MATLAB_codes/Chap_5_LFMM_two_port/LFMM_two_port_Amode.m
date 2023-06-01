%% LFMM two port analysis

%case 1. free space - finite PC waveguide - free space 
%case 2. free space - half-infinite PC waveguide 
%case 3. half-infinite PC waveguide - free space

%% STEP 1
clear all;
close all;
clc;

addpath([pwd '\PRCWA_COM']);
addpath([pwd '\FIELD_VISUAL']);
addpath([pwd '\STRUCTURE']);

% length unit

global nm; global nano;     % nano
global um; global micro;    % micro
global mm;                  % mili
global lambda;              % wavelength

global k0;                                  % wavenumber
global c0; global w0;
global eps0; global mu0;

% zero-thickness buffer refractive index, permittivity, permeability
global n0; global epr0; global mur0; 
% refractive index, permittivity, permeability in free space I
global ni; global epri; global muri;
% refractive index, permittivity, permeability in free space III
global nf; global eprf; global murf;

% x-directional supercell period, y-directional supercell period 
global aTx; global aTy; global aTz; global nx; global ny; global nz;
global bTx; global bTy; global bTz; global Tx; global Ty; global Tz;
global NBx; global NBy; global NBz; global num_hx; global num_hy; global k0;
global kx_vc; global ky_vc; global kz_vc;

% input output free space
global kix; global kiy; global kiz; global kfz;
global kx_ref; global ky_ref; global kz_ref;
global kx_tra; global ky_tra; global kz_tra;

PRCWA_basic;        % 3D structure
PRCWA_Gen_K;        % zero-thickness buffer
PRCWA_Gen_inout_Ka;  %          port1-port2 input-output region     

% structure modeling
% The example structure is the photonic crystal waveguide structure
% approxiated by the stepwise multi-blocks. The structure modeling is
% implemented in the MATLAB code line 59, PRCWA_Gen_PCWG.m
% In PRCWA_Gen_PCWG.m, the Fourier series of the target structure is done
% by Grating_gen_Y_branch.m as Grating_Gen_PCWG.m. Open PRCWA_Gen_PCWG.m
% and see the annotation

PRCWA_Gen_PCWG; % photonic crystal waveguide

%% STEP2 Block S-matrix computation of single super-block structures
% Bloch eigenmode analysis is performed by AmodeBlochAnalysis.m. Let us see
% its inside in code

%%% A mode Bloch Mode Analysis
% Fourier coefficients of A Bloch eigenmodes

ApEFx=zeros(NBx,NBy,NBz,2*L);
ApEFy=zeros(NBx,NBy,NBz,2*L);
ApEFz=zeros(NBx,NBy,NBz,2*L);

ApHFx=zeros(NBx,NBy,NBz,2*L);
ApHFy=zeros(NBx,NBy,NBz,2*L);
ApHFz=zeros(NBx,NBy,NBz,2*L);

AnEFx=zeros(NBx,NBy,NBz,2*L);
AnEFy=zeros(NBx,NBy,NBz,2*L);
AnEFz=zeros(NBx,NBy,NBz,2*L);

AnHFx=zeros(NBx,NBy,NBz,2*L);
AnHFy=zeros(NBx,NBy,NBz,2*L);
AnHFz=zeros(NBx,NBy,NBz,2*L);

% Positive & Negative eigenvalue of a modes
ap_evalue=zeros(1,2*L);
am_evalue=zeros(1,2*L);

AmodeBlochAnalysis;

ApEFx=pEFx;
ApEFy=pEFy;
ApEFz=pEFz;

ApHFx=pHFx;
ApHFy=pHFy;
ApHFz=pHFz;

AnEFx=nEFx;
AnEFy=nEFy;
AnEFz=nEFz;

AnHFx=nHFx;
AnHFy=nHFy;
AnHFz=nHFz;

ap_evalue=pe_value;
am_evalue=me_value;
ap_evector=pe_vector;
am_evector=me_vector;

AModePower_Flow;

%% STEP 3 S-matrix computation

% Eqs. (5.1.7a)~(5.1.7d) can be calculated efficiently with the use of
% bultin DFT function of MATLAB. The field calculation of Bloch eigenmodes
% is implemented in Amode_visualization_fft.m

% A mode [2-port Block]  
% Layer S-matrix of Amode
AR11=zeros(2*L,2*L); %  port1->port1
AT12=zeros(2*L,2*L); %  port1->port2
AT21=zeros(2*L,2*L); %  port2->port1
AR22=zeros(2*L,2*L); %  port2->port2
ACa=zeros(4*L,2*L); %   excitation of port1
ACb=zeros(4*L,2*L); %   excitation of port2

% Left grating S-matrix
ALR11=zeros(2*L,2*L); % port1->port1
ALT12=zeros(2*L,2*L); % port1->port2
ALT21=zeros(2*L,2*L); % port2->port1
ALR22=zeros(2*L,2*L); % port2->port2

% Right grating S-matrix
ARR11=zeros(2*L,2*L); %  port1->port1
ART12=zeros(2*L,2*L); % port1->port2
ART21=zeros(2*L,2*L); % port2->port1
ARR22=zeros(2*L,2*L); % port2->port2

zm=0;
zp=aTz;
zc=0;

SMM_2port_Amode; 

% A mode eigenmode calculation
Amode_visualization_fft;

%% STEP4 Half-infinite block interconnection & Field visualization

LFMM_Smat_Amode_case1;          % case1
%LFMM_Smat_Amode_case2;         % case2
%LFMM_Smat_Amode_case3;         % case3

%% Data save
save 2portDATA_Amode_all;
save 2portDATA_Amode aPG_Ex_xz aPG_Ey_xz aPG_Ez_xz...
     aPG_Hx_xz aPG_Hy_xz aPG_Hz_xz...
     aNG_Ex_xz aNG_Ey_xz aNG_Ez_xz...
     aNG_Hx_xz aNG_Hy_xz aNG_Hz_xz...
     ALR11 ALT12 ALR22 ALT21...
     ARR11 ART12 ARR22 ART21...
     AR11  AT12  AT21  AR22 ACa ACb...
     ap_evalue am_evalue ap_evector am_evector...
     aPw aNw aTx aTz SAR11 SAT12 SAT21 SAR22 SACa SACb;
 