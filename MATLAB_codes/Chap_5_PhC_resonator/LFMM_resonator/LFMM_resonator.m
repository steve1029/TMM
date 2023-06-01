%% LFMM with the staircase approximation
% 2010 01 04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LFMM with the extended S-matrix method
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
tic


addpath([pwd '\PRCWA_COM']);
addpath([pwd '\FIELD_VISUAL']);
addpath([pwd '\STRUCTURE']);


%% part 1-2

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


%% part 2

PRCWA_basic;        % 3D structure
PRCWA_Gen_K;        % zero-thickness buffer
PRCWA_Gen_inout_Ka;  %          port1-port2 input-output region
%PRCWA_Gen_inout_Kb;  %          port3-port4 input-output region

% Positive & Negative eigenvalue of a modes
ap_evalue=zeros(1,2*L);
am_evalue=zeros(1,2*L);

% Positive & Negative eigenvalue of b modes
% bp_evalue=zeros(1,2*L);
% bm_evalue=zeros(1,2*L);


%-----------------------------------------------------
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

% Fourier coefficients of B Bloch eigenmodes

% BpEFx=zeros(NBz,NBy,NBx,2*L);
% BpEFy=zeros(NBz,NBy,NBx,2*L);
% BpEFz=zeros(NBz,NBy,NBx,2*L);
% 
% BpHFx=zeros(NBz,NBy,NBx,2*L);.
% BpHFy=zeros(NBz,NBy,NBx,2*L);
% BpHFz=zeros(NBz,NBy,NBx,2*L);
% 
% BnEFx=zeros(NBz,NBy,NBx,2*L);
% BnEFy=zeros(NBz,NBy,NBx,2*L);
% BnEFz=zeros(NBz,NBy,NBx,2*L);
% 
% BnHFx=zeros(NBz,NBy,NBx,2*L);
% BnHFy=zeros(NBz,NBy,NBx,2*L);
% BnHFz=zeros(NBz,NBy,NBx,2*L);


%% part 3: structure modeling

% PRCWA_Gen_PCWG; % photonic crystal waveguide
%PRCWA_Gen_PhC_waveguide
PRCWA_Gen_PhC_resonator

%% part 4-1 Layer S-matrix computation of layers


%%% A mode Bloch Mode Analysis
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


%%% B mode Bloch Mode Analysis
% BmodeBlochAnalysis;
% for k=1:NBx
%     for l=1:NBy
%         for m=1:NBz
%             BpEFx(m,l,k,:)=pEFz(NBx-k+1,l,m,:);     %%%% very important!!!
%             BpEFy(m,l,k,:)=pEFy(NBx-k+1,l,m,:);
%             BpEFz(m,l,k,:)=-pEFx(NBx-k+1,l,m,:);
% 
%             BpHFx(m,l,k,:)=pHFz(NBx-k+1,l,m,:);
%             BpHFy(m,l,k,:)=pHFy(NBx-k+1,l,m,:);
%             BpHFz(m,l,k,:)=-pHFx(NBx-k+1,l,m,:);
% 
%             BnEFx(m,l,k,:)=nEFz(NBx-k+1,l,m,:);
%             BnEFy(m,l,k,:)=nEFy(NBx-k+1,l,m,:);
%             BnEFz(m,l,k,:)=-nEFx(NBx-k+1,l,m,:);
% 
%             BnHFx(m,l,k,:)=nHFz(NBx-k+1,l,m,:);
%             BnHFy(m,l,k,:)=nHFy(NBx-k+1,l,m,:);
%             BnHFz(m,l,k,:)=-nHFx(NBx-k+1,l,m,:);
%         end
%     end
% end
% 
% bp_evalue=pe_value;
% bm_evalue=me_value;
% bp_evector=pe_vector;
% bm_evector=me_vector;
% BModePower_Flow;

%% S-matrix computation

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

% % B mode [2-port Block]
% % Layer S-matrix of Bmode
% BR33=zeros(2*L,2*L); %  port3->port3
% BT34=zeros(2*L,2*L); %  port3->port4
% BT43=zeros(2*L,2*L); %  port4->port3
% BR44=zeros(2*L,2*L); %  port4->port4
% BCa=zeros(4*L,2*L); %   excitation of port3
% BCb=zeros(4*L,2*L); %   excitation of port4
% 
% % Lower Boundary S-matrix
% BLR33=zeros(2*L,2*L); % port3->port3
% BLT34=zeros(2*L,2*L); % port3->port4
% BLT43=zeros(2*L,2*L); % port4->port3
% BLR44=zeros(2*L,2*L); % port4->port4
% 
% % Upper Boundary S-matrix
% BRR33=zeros(2*L,2*L); % port3->port3
% BRT34=zeros(2*L,2*L); % port3->port4
% BRT43=zeros(2*L,2*L); % port4->port3
% BRR44=zeros(2*L,2*L); % port4->port4
% 
% xm=0;
% xp=bTx;
% xc=0;
% 
% % SMM_2port_Bmode;


%% Field visualization
%case 1. homogeneous space          - finite PC waveguide - homogeneous space
%case 2. homogeneous space          - semi-infinite PC waveguide
%case 3. semi-infinite PC waveguide - homogeneous space

% A mode
Amode_visualization_fft;
toc
disp('done~!');
%return;
%LFMM_Smat_Amode_case1;        % case1
%LFMM_Smat_Amode_case2;        % case2
%LFMM_Smat_Amode_case3;         % case3


% % B mode
% Bmode_visualization_fft;
% 
% % LFMM_Smat_Bmode_case1;         % case1
% % LFMM_Smat_Bmode_case2;         % case2
% % LFMM_Smat_Bmode_case3;         % case3

save 2portDATA_resonator;

% save 2portDATA aPG_Ex_xz aPG_Ey_xz aPG_Ez_xz...
%     aPG_Hx_xz aPG_Hy_xz aPG_Hz_xz...
%     aNG_Ex_xz aNG_Ey_xz aNG_Ez_xz...
%     aNG_Hx_xz aNG_Hy_xz aNG_Hz_xz...
%     bPG_Ex_xz bPG_Ey_xz bPG_Ez_xz...
%     bPG_Hx_xz bPG_Hy_xz bPG_Hz_xz...
%     bNG_Ex_xz bNG_Ey_xz bNG_Ez_xz...
%     bNG_Hx_xz bNG_Hy_xz bNG_Hz_xz...
%     ALR11 ALT12 ALR22 ALT21...
%     ARR11 ART12 ARR22 ART21...
%     AR11  AT12  AT21  AR22 ACa ACb...
%     BLR33 BLT34 BLR44 BLT43...
%     BRR33 BRT34 BRR44 BRT43...
%     BR33  BT34 BT43  BR44 BCa  BCb...
%     ap_evalue am_evalue bp_evalue bm_evalue...
%     ap_evector am_evector bp_evector bm_evector...
%     aPw aNw bPw bNw aTx aTz bTx bTz...
%     SAR11 SAT12 SAT21 SAR22 SACa SACb ...
%     SBR33 SBT34 SBT43 SBR44 SBCa SBCb;








