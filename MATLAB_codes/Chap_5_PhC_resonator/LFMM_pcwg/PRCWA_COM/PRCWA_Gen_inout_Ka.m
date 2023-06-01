% PRCWA_Gen_inout_Ka

% left & right direction wavevector components
aki=k0*ni;                      % left input region wavenumber
akf=k0*nf;                      % right output region wavenumber

akix=aki*sin(theta)*cos(phi);   % incident wavevector kx
akiy=aki*sin(theta)*sin(phi);   % incident wavevector ky
akiz=aki*cos(theta);            % incident wavevector kz
akfz=sqrt(akf^2-akix^2-akiy^2); % output region wavevector kz

% region I : incident & reflection 
akx_ref=zeros(NBx,1);    % supported mode number in x direction
aky_ref=zeros(NBy,1);    % supported mode number in y direction
akz_ref=zeros(NBx,NBy);  % supported mode number in z direction
R=zeros(NBx*NBy,1);      

for k=1:NBx
for l=1:NBy  
	akx_ref(k)=aki*sin(theta)*cos(phi)+(k-(nx+1))*(2*pi/Tx); 
	aky_ref(l)=aki*sin(theta)*sin(phi)+(l-(ny+1))*(2*pi/Ty); 
	akz_ref(k,l)=k0*(ni^2-(akx_ref(k)/k0)^2-(aky_ref(l)/k0)^2)^0.5;	
end;
end;

% region II : transmission
akx_tra=zeros(NBx,1);  	 % supported mode number in x direction
aky_tra=zeros(NBy,1);  	 % supported mode number in y direction
akz_tra=zeros(NBx,NBy);  % supported mode number in z direction

T=zeros(NBx*NBy,1);      
for k=1:NBx
for l=1:NBy
   akx_tra(k)=akx_ref(k);
   aky_tra(l)=aky_ref(l);
   akz_tra(k,l)=k0*(nf^2-(akx_ref(k)/k0)^2-(aky_ref(l)/k0)^2)^0.5;
end;    % for l
end;	% for k
