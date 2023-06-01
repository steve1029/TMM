% odd sampling fftshift

function y=odd_fftshift(x,nx,ny)

tx=zeros(nx,ny);
ttx=zeros(nx,ny);

tx=x_fftshift(x,nx,ny);
ttx=z_fftshift(tx,nx,ny);
y=ttx;

%for k=1:(ny-1)/2
%   tx(:,k+1)=x(:,k+(ny+1)/2);
%    tx(:,k+(ny+1)/2)=x(:,k);
%end;
%    tx(:,1)=x(:,(ny+1)/2);
    
%for k=1:(nx-1)/2
%    ttx(k+1,:)=tx(k+(nx+1)/2,:);
%    ttx(k+(nx+1)/2,:)=tx(k,:);
%end;
%    ttx(1,:)=tx((nx+1)/2,:);
    
%y=ttx;