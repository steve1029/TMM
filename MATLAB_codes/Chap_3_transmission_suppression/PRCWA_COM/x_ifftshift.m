% x-direction ifftshift
function y=x_ifftshift(x,nx,ny)

tx=zeros(nx,ny);

for k=1:(nx-1)/2
    tx(k,:)=x(k+(nx+1)/2,:);
    tx(k+(nx-1)/2,:)=x(k,:);
end;
    tx(nx,:)=x((nx+1)/2,:);
    y=tx;