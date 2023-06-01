% rect_2D_mesh.m

function y=rect_2D_mesh(k,l,ep,Tx,Ty,x1,x2,y1,y2)

dx=abs(x1-x2);
dy=abs(y1-y2);


f1=( (dx/Tx)*sinc(dx*k/Tx).*exp(-j*pi*k/Tx*(x1+x2)) );
f2=( (dy/Ty)*sinc(dy*l/Ty).*exp(-j*pi*l/Ty*(y1+y2)) );
y=f1.*f2;
