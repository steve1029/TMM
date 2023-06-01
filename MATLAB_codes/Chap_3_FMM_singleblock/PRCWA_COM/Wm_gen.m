function Wm_z=Wm_gen(mW,mevalue,mcnt,L,z)

global c0; global w0;
global eps0; global mu0;

Wm_z=zeros(2*L,mcnt);

for k=1:mcnt

Wm_z(:,k)=mW(1:2*L,k)*exp(mevalue(k)*z);

end;