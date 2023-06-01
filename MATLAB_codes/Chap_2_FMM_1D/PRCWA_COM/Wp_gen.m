function Wp_z=Wp_gen(pW,pevalue,pcnt,L,z)

global c0; global w0;
global eps0; global mu0;

Wp_z=zeros(2*L,pcnt);

for k=1:pcnt

Wp_z(:,k)=pW(1:2*L,k)*exp(pevalue(k)*z);

end;

