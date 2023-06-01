% PRCWA_Gen_inout_K
% left & right wavevector components
ki=k0*ni;                   % left input   region wavenumber
kf=k0*nf;                   % right output region wavenumber

kix=ki*sin(theta)*cos(phi); % incident wavevector kx
kiy=ki*sin(theta)*sin(phi); % incident wavevector ky
kiz=(ki^2-kix^2-kiy^2)^0.5; % incident wavevector kz
kfz=(kf^2-kix^2-kiy^2)^0.5; % output region wavevector kz

% region I : incident & reflection 
kx_ref=zeros(NBx,1);    % supported mode number in x direction: kxm
ky_ref=zeros(NBy,1);    % supported mode number in y direction: kxn
kz_ref=zeros(NBx,NBy);  % supported mode number in z direction: kzmn
R=zeros(NBx*NBy,1);     

for k=1:NBx
    for l=1:NBy  
        kx_ref(k)=kix+(k-(nx+1))*(2*pi/Tx);          
        ky_ref(l)=kiy+(l-(ny+1))*(2*pi/Ty); 
        kz_ref(k,l)=k0*(ni^2-(kx_ref(k)/k0)^2-(ky_ref(l)/k0)^2)^0.5;	
    end;
end;

L=NBx*NBy;


KI=zeros(2*L,2*L); % KI
for k=1:NBx
   for l=1:NBy
od_ind1=(k-1)*NBy+l;      
K11(od_ind1,od_ind1)=(kx_ref(k)*ky_ref(l))/kz_ref(k,l);
K12(od_ind1,od_ind1)=(kz_ref(k,l)^2+kx_ref(k)^2)/kz_ref(k,l);
K21(od_ind1,od_ind1)=-(ky_ref(l)^2+kz_ref(k,l)^2)/kz_ref(k,l);
K22(od_ind1,od_ind1)=-( ky_ref(l)*kx_ref(k) )/kz_ref(k,l);

   end;
end;

KI(1:L, 1:L)=K11;
KI(1:L, L+1:2*L)=K12;
KI(L+1:2*L, 1:L)=K21;
KI(L+1:2*L,L+1:2*L )=K22;


% region II : transmission
kx_tra=zeros(NBx,1);  	 % supported mode number in x direction
ky_tra=zeros(NBy,1);  	 % supported mode number in y direction
kz_tra=zeros(NBx,NBy);   % supported mode number in z direction
T=zeros(NBx*NBy,1); 

for k=1:NBx
    for l=1:NBy
        kx_tra(k)=kx_ref(k);
        ky_tra(l)=ky_ref(l);
        kz_tra(k,l)=k0*(nf^2-(kx_ref(k)/k0)^2-(ky_ref(l)/k0)^2)^0.5;
    end; % for l
end;	% for k

KII=zeros(2*L,2*L); % KII
for k=1:NBx
   for l=1:NBy
od_ind1=(k-1)*NBy+l;      
K11(od_ind1,od_ind1)=(kx_tra(k)*ky_tra(l))/kz_tra(k,l);
K12(od_ind1,od_ind1)=(kz_tra(k,l)^2+kx_tra(k)^2)/kz_tra(k,l);
K21(od_ind1,od_ind1)=-(ky_tra(l)^2+kz_tra(k,l)^2)/kz_tra(k,l);
K22(od_ind1,od_ind1)=-(ky_tra(l)*kx_tra(k) )/kz_tra(k,l);

   end;
end;

KII(1:L, 1:L)=K11;
KII(1:L, L+1:2*L)=K12;
KII(L+1:2*L, 1:L)=K21;
KII(L+1:2*L,L+1:2*L )=K22;


Wh=Iden;
VhI=KI/(w0*mu0);
VhII=KII/(w0*mu0);
