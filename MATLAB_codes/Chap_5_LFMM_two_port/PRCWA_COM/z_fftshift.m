% z-direction fftshift

function y=z_fftshift(x,nx,ny)

tx=zeros(nx,ny);

for k=1:(ny-1)/2
    tx(:,k+1)=x(:,k+(ny+1)/2);
    tx(:,k+(ny+1)/2)=x(:,k);
end;

    tx(:,1)=x(:,(ny+1)/2);
    
 y=tx;

