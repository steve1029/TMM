% rect_1D_mesh


function y=rect_1D_mesh(k,ep,d,l1,l2)
dp=abs(l1-l2);
y=(ep*dp/d)*sinc(dp*k/d).*exp(-j*pi*k/d*(l1+l2));