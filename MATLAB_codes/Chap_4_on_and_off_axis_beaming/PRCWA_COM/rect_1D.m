% rect_1D


function y=rect_1D(k,ep,d,l1,l2)
dp=abs(l1-l2);
y=(ep*dp/d)*sinc(dp*k/d)*exp(-j*pi*k/d*(l1+l2));