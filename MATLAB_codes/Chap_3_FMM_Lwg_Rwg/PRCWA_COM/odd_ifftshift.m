% odd sampling fftshift

function y=odd_ifftshift(x,nx,ny)

tx=zeros(nx,ny);
ttx=zeros(nx,ny);

tx=x_ifftshift(x,nx,ny);
ttx=z_ifftshift(tx,nx,ny);
y=ttx;
