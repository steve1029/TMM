function Vp_z=Vp_gen(pW,pevalue,pcnt,L,z)

global c0; global w0;
global eps0; global mu0;

Vp_z=zeros(2*L,pcnt);

for k=1:pcnt
      
Vp_z(:,k)=j*(eps0/mu0)^0.5*pW(2*L+1:4*L,k)*exp(pevalue(k)*z);

end;

