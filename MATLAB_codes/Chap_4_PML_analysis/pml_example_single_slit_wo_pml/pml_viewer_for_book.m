xx=(-Tx/2:Tx*0.001:Tx/2);
Gr_str=zeros(1,length(xx));

for k=-2*nx:2*nx
   Gr_str=Gr_str+Epsr_PML(k+NBx)*exp(j*(k*xx*2*pi/Tx));
end
figure(31);set(gca,'fontsize',16);set(gca,'fontname','times new roman');
plot(real(Epsr_PML));
figure(32);set(gca,'fontsize',16);set(gca,'fontname','times new roman');
plot(xx/Tx,real(Gr_str),'r','linewidth',2);hold on;set(gca,'fontname','times new roman');
xlabel('x (in unit of T_x)');set(gca,'fontname','times new roman');
ylabel('\epsilon');set(gca,'fontname','times new roman');
plot(xx/Tx,imag(Gr_str),':b','linewidth',2);
