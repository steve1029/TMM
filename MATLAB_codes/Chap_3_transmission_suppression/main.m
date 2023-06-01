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

addpath([pwd '\PRCWA_COM']);
addpath([pwd '\FIELD_VISUAL']);
addpath([pwd '\STRUCTURE']);


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

nm=1e-9;
lambda=514.5*nm;

Tx_set=(150:5:950);Tx_len=length(Tx_set);
result_trans=zeros(2*(50)+1,Tx_len);
result_refle=zeros(2*(50)+1,Tx_len);
result_cetHy=zeros(1,Tx_len);

t1=clock;
t2=clock;
for Tx_idx = 1 : Tx_len
    Tx=(Tx_set(Tx_idx)+1e-2) * nm;Ty=Tx;
    accum_time=etime(t2,t1);
    if Tx_idx ~= 1
        remaining=(Tx_len-Tx_idx+1)*accum_time/(Tx_idx-1);
    else
        remaining=0;
    end
    disp([num2str(Tx_idx) '/' num2str(Tx_len)  ' ela : ' num2str(accum_time) ' rem : ' num2str(remaining)]);


%% part 2

PRCWA_basic;                                %   3D structure
PRCWA_Gen_K;                                %   zero-thickness buffer

%% part 3 : structure modeling
%PRCWA_Gen_PML2D;                            %   PML build
%PRCWA_Gen_SPP_beaming;                     %   SPP beaming
% PRCWA_Gen_Y_branch;                         %   SPP Y branch
PRCWA_Gen_diagonal_TriangleGrating2;         %   EOT structure

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


alpha_tm=1;
beta_tm=1;
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
 
%% part 4-2 S-matrix method of boundary matching
 
 Bdr_Smat_case1;      % 1. homogeneous space          - grating - homogeneous space    
 %Bdr_Smat_case2;       % 2. homogeneous space          - grating - inhomogeneous waveguide
 %Bdr_Smat_case3;      % 3. inhomogeneous waveguide    - grating - homogeneous space
%  Bdr_Smat_case4;      % 4. inhomogeneous waveguide    - grating - inhomogeneous waveguide

%% part 5 Field visualization & data save

Field_visualization_3D_xz_case1_Lfree_Rfree_leftright;
% Field_visualization_3D_xz_case1_Lfree_Rfree_rightleft;

% Field_visualization_3D_yz_case1_Lfree_Rfree_leftright;
% Field_visualization_3D_yz_case1_Lfree_Rfree_rightleft;

%% part 6 save the results

result_refle(:,Tx_idx)=DEt1;
result_trans(:,Tx_idx)=DEt3;
t2=clock;

end

%%
min_val=min(sum(result_trans));
max_val=max(sum(result_refle));

figure(11);set(gca,'fontsize',16);set(gca,'fontname','times new roman');box on;
semilogy(Tx_set,sum(result_refle),':r','linewidth',2);hold on;
semilogy(Tx_set,sum(result_trans),'-b','linewidth',2);
axis([Tx_set(1) Tx_set(end) min_val*1.0 max_val*1.1]);set(gca,'fontname','times new roman');
xlabel('Grating period (nm)');set(gca,'fontname','times new roman');
legend('Reflection','Transmission');set(gca,'fontname','times new roman');

lambda_spp=307;
line(1*lambda_spp*ones(1,10),linspace(min_val,max_val,10),'linestyle','--','color','k');
line(2*lambda_spp*ones(1,10),linspace(min_val,max_val,10),'linestyle','--','color','k');
line(3*lambda_spp*ones(1,10),linspace(min_val,max_val,10),'linestyle','--','color','k');


