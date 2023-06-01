% PRCWA_Gen_inout_K
% left & right wavevector components
ki=k0*ni;                   % left input   region wavenumber
kf=k0*nf;                   % right output region wavenumber

kix=kx0; % incident wavevector kx
kiy=ky0; % incident wavevector ky
kiz=(ki^2-kix^2-kiy^2)^0.5; % incident wavevector kz
kfz=(kf^2-kix^2-kiy^2)^0.5; % output region wavevector kz

% region I : incident & reflection 

        kx_ref=kix;          
        ky_ref=kiy; 
        kz_ref=k0*(ni^2-(kx_ref/k0)^2-(ky_ref/k0)^2)^0.5;	
  
L=1;

KI=zeros(2*L,2*L); % KI
     
K11=(kx_ref*ky_ref)/kz_ref;
K12=(kz_ref^2+kx_ref^2)/kz_ref;
K21=-(ky_ref^2+kz_ref^2)/kz_ref;
K22=-( ky_ref*kx_ref)/kz_ref;

KI(1:L, 1:L)=K11;
KI(1:L, L+1:2*L)=K12;
KI(L+1:2*L, 1:L)=K21;
KI(L+1:2*L,L+1:2*L )=K22;


% region II : transmission

        kx_tra=kx_ref;
        ky_tra=ky_ref;
        kz_tra=k0*(nf^2-(kx_ref/k0)^2-(ky_ref/k0)^2)^0.5;


KII=zeros(2*L,2*L); % KII
      
K11=(kx_tra*ky_tra)/kz_tra;
K12=(kz_tra^2+kx_tra^2)/kz_tra;
K21=-(ky_tra^2+kz_tra^2)/kz_tra;
K22=-(ky_tra*kx_tra)/kz_tra;

KII(1:L, 1:L)=K11;
KII(1:L, L+1:2*L)=K12;
KII(L+1:2*L, 1:L)=K21;
KII(L+1:2*L,L+1:2*L )=K22;


Wh=Iden;
VhI=KI/(w0*mu0);
VhII=KII/(w0*mu0);
