% x-direction fftshift

function y=x_fftshift(x,nx,ny)

tx=zeros(nx,ny);

for k=1:(nx-1)/2
    tx(k+1,:)=x(k+(nx+1)/2,:);
    tx(k+(nx+1)/2,:)=x(k,:);
end;
    tx(1,:)=x((nx+1)/2,:);
    
    y=tx;