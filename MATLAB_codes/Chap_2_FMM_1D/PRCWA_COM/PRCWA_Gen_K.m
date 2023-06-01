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






