function Vm_z=sVm_gen(mV,mevalue,mcnt,L,z)

global c0; global w0;
global eps0; global mu0;

Vm_z=zeros(2*L,mcnt);

for k=1:mcnt

Vm_z(:,k)=j*(eps0/mu0)^0.5*mV(1:2*L,k)*exp(mevalue(k)*z);

end;

