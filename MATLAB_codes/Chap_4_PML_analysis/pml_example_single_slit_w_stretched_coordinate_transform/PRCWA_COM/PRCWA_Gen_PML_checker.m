xx=1*(-Tx/2:Tx*0.01:Tx/2)';
Gr_str_E00=zeros(length(xx),Nlay);
Gr_str_Exx=zeros(length(xx),Nlay);
Gr_str_Eyy=zeros(length(xx),Nlay);
Gr_str_Ezz=zeros(length(xx),Nlay);
Gr_str_Axx=zeros(length(xx),Nlay);
Gr_str_Ayy=zeros(length(xx),Nlay);
Gr_str_Azz=zeros(length(xx),Nlay);

% for laynt=1:Nlay
laynt=3;
    for k=-2*nx:2*nx
%         Gr_str_Exx(:,laynt)=Gr_str_Exx(:,laynt)+eps_xx(k+NBx,laynt)*exp(j*(k*xx*2*pi/Tx));
%         Gr_str_Eyy(:,laynt)=Gr_str_Eyy(:,laynt)+eps_yy(k+NBx,laynt)*exp(j*(k*xx*2*pi/Tx));
%         Gr_str_Ezz(:,laynt)=Gr_str_Ezz(:,laynt)+eps_zz(k+NBx,laynt)*exp(j*(k*xx*2*pi/Tx));
%         Gr_str_Axx(:,laynt)=Gr_str_Axx(:,laynt)+aps_xx(k+NBx,laynt)*exp(j*(k*xx*2*pi/Tx));
%         Gr_str_Ayy(:,laynt)=Gr_str_Ayy(:,laynt)+aps_yy(k+NBx,laynt)*exp(j*(k*xx*2*pi/Tx));
%         Gr_str_Azz(:,laynt)=Gr_str_Azz(:,laynt)+aps_zz(k+NBx,laynt)*exp(j*(k*xx*2*pi/Tx));
        
        Gr_str_E00(:,laynt)=Gr_str_E00(:,laynt)+rect_region0(k+NBx,1)*exp(j*(k*xx*2*pi/Tx));

        Gr_str_Exx(:,laynt)=Gr_str_Exx(:,laynt)+Epsr_PML_xx(k+NBx,1)*exp(j*(k*xx*2*pi/Tx));
        Gr_str_Eyy(:,laynt)=Gr_str_Eyy(:,laynt)+Epsr_PML_yy(k+NBx,1)*exp(j*(k*xx*2*pi/Tx));
        Gr_str_Ezz(:,laynt)=Gr_str_Ezz(:,laynt)+Epsr_PML_zz(k+NBx,1)*exp(j*(k*xx*2*pi/Tx));
        Gr_str_Axx(:,laynt)=Gr_str_Axx(:,laynt)+Apsr_PML_xx(k+NBx,1)*exp(j*(k*xx*2*pi/Tx));
        Gr_str_Ayy(:,laynt)=Gr_str_Ayy(:,laynt)+Apsr_PML_yy(k+NBx,1)*exp(j*(k*xx*2*pi/Tx));
        Gr_str_Azz(:,laynt)=Gr_str_Azz(:,laynt)+Apsr_PML_zz(k+NBx,1)*exp(j*(k*xx*2*pi/Tx));
    end
% end

%%
figure(71);
subplot(1,2,1);plot(real(Gr_str_E00));
subplot(1,2,2);plot(imag(Gr_str_E00));


%%
figure(50+laynt);
subplot(2,3,1);plot(real(Gr_str_Exx));
subplot(2,3,2);plot(real(Gr_str_Eyy));
subplot(2,3,3);plot(real(Gr_str_Ezz));
subplot(2,3,4);plot(real(Gr_str_Axx));
subplot(2,3,5);plot(real(Gr_str_Ayy));
subplot(2,3,6);plot(real(Gr_str_Azz));

%%
figure(60+laynt);
subplot(2,3,1);plot(imag(Gr_str_Exx));
subplot(2,3,2);plot(imag(Gr_str_Eyy));
subplot(2,3,3);plot(imag(Gr_str_Ezz));
subplot(2,3,4);plot(imag(Gr_str_Axx));
subplot(2,3,5);plot(imag(Gr_str_Ayy));
subplot(2,3,6);plot(imag(Gr_str_Azz));

