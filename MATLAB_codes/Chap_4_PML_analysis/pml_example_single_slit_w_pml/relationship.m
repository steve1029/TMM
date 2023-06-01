em=-10.1798 + 0.8259*i;
ea=1;
ed=1.72^2;

nspp_a=real(sqrt(em*ea/(em+ea)));
nspp_d=real(sqrt(em*ed/(em+ed)));
% nspp=(nspp_a+nspp_d)/2;
nspp=1.4;

wavelen=532;

periods=linspace(250,550,101);

thetas=asind(nspp-wavelen./periods);

figure(3);set(gca,'fontsize',16);set(gca,'fontname','times new roman');
plot(periods,thetas,'k');grid on;
axis([periods(1) periods(end) -30 30]);set(gca,'fontname','times new roman');
xlabel('Grating period (nm)');
ylabel('Diffraction angle (deg.)');