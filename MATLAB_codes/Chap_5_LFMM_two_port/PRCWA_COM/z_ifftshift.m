% z-direction ifftshift
function y=z_ifftshift(x,nx,ny)

tx=zeros(nx,ny);

for k=1:(ny-1)/2
    tx(:,k)=x(:,k+(ny+1)/2);
    tx(:,k+(ny-1)/2)=x(:,k);
end;

    tx(:,ny)=x(:,(ny+1)/2);
    
 y=tx;

