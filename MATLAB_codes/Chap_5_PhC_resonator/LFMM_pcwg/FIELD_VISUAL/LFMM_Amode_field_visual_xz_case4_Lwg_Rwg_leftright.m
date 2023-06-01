%LFMM_Amode_field_visual_xz_case4_Lwg_Rwg_leftright


D=zeros(2*L,1);
inc_mode=1;
D(inc_mode)=1;

%region I
CR1=RRa*D;
Lay_cof1=zeros(4*L,1);
Lay_cof1(1:2*L,1)=D;
Lay_cof1(2*L+1:4*L,1)=CR1;

% region III
CT3=TTa*D;
Lay_cof3=zeros(4*L,1);
Lay_cof3(1:2*L,1)=CT3;

% grating : region II
SNlay_resonator=1;

clear C_temp;
clear C_p;
clear C_n;

for laynt=1:SNlay_resonator
    C_temp(:,laynt)=Ca(:,:,laynt)*D;
    C_p(:,laynt)=C_temp(1:2*L,laynt);          % Positive coupling coefficients
    C_n(:,laynt)=C_temp(2*L+1:4*L,laynt);      % Negative coupling coefficients
end;



%% field visualization - region I

SNlay_wg1=1;

[xlen zlen]=size(aNG_Ey_xz_wg(:,:,1));  
aGEy_xz1=zeros(xlen,zlen,SNlay_wg1);
aGEx_xz1=zeros(xlen,zlen,SNlay_wg1);
aGEz_xz1=zeros(xlen,zlen,SNlay_wg1);
aGHy_xz1=zeros(xlen,zlen,SNlay_wg1);
aGHx_xz1=zeros(xlen,zlen,SNlay_wg1);
aGHz_xz1=zeros(xlen,zlen,SNlay_wg1);

for laynt=1:SNlay_wg1
    for pbm_ind=1:mcnt
        
        aGEy_xz1(:,:,laynt)=aGEy_xz1(:,:,laynt)+CR1(pbm_ind)*aNG_Ey_xz_wg(:,:,pbm_ind)*exp(-am_evalue_wg(pbm_ind)*aTz_wg*(SNlay_wg1-laynt)); %Sx
        aGEx_xz1(:,:,laynt)=aGEx_xz1(:,:,laynt)+CR1(pbm_ind)*aNG_Ex_xz_wg(:,:,pbm_ind)*exp(-am_evalue_wg(pbm_ind)*aTz_wg*(SNlay_wg1-laynt)); %Sx
        aGEz_xz1(:,:,laynt)=aGEz_xz1(:,:,laynt)+CR1(pbm_ind)*aNG_Ez_xz_wg(:,:,pbm_ind)*exp(-am_evalue_wg(pbm_ind)*aTz_wg*(SNlay_wg1-laynt));
        
        aGHy_xz1(:,:,laynt)=aGHy_xz1(:,:,laynt)+CR1(pbm_ind)*aNG_Hy_xz_wg(:,:,pbm_ind)*exp(-am_evalue_wg(pbm_ind)*aTz_wg*(SNlay_wg1-laynt)); %Sx
        aGHx_xz1(:,:,laynt)=aGHx_xz1(:,:,laynt)+CR1(pbm_ind)*aNG_Hx_xz_wg(:,:,pbm_ind)*exp(-am_evalue_wg(pbm_ind)*aTz_wg*(SNlay_wg1-laynt)); %Sx
        aGHz_xz1(:,:,laynt)=aGHz_xz1(:,:,laynt)+CR1(pbm_ind)*aNG_Hz_xz_wg(:,:,pbm_ind)*exp(-am_evalue_wg(pbm_ind)*aTz_wg*(SNlay_wg1-laynt));
        
    end;
    % incidence mode
    aGEy_xz1(:,:,laynt)=aGEy_xz1(:,:,laynt)+aPG_Ey_xz_wg(:,:,inc_mode)*exp(-ap_evalue_wg(inc_mode)*aTz_wg*(SNlay_wg1-laynt+1)); %Sx
    aGEx_xz1(:,:,laynt)=aGEx_xz1(:,:,laynt)+aPG_Ex_xz_wg(:,:,inc_mode)*exp(-ap_evalue_wg(inc_mode)*aTz_wg*(SNlay_wg1-laynt+1)); %Sx
    aGEz_xz1(:,:,laynt)=aGEz_xz1(:,:,laynt)+aPG_Ez_xz_wg(:,:,inc_mode)*exp(-ap_evalue_wg(inc_mode)*aTz_wg*(SNlay_wg1-laynt+1));
    
    aGHy_xz1(:,:,laynt)=aGHy_xz1(:,:,laynt)+aPG_Hy_xz_wg(:,:,inc_mode)*exp(-ap_evalue_wg(inc_mode)*aTz_wg*(SNlay_wg1-laynt+1)); %Sx
    aGHx_xz1(:,:,laynt)=aGHx_xz1(:,:,laynt)+aPG_Hx_xz_wg(:,:,inc_mode)*exp(-ap_evalue_wg(inc_mode)*aTz_wg*(SNlay_wg1-laynt+1)); %Sx
    aGHz_xz1(:,:,laynt)=aGHz_xz1(:,:,laynt)+aPG_Hz_xz_wg(:,:,inc_mode)*exp(-ap_evalue_wg(inc_mode)*aTz_wg*(SNlay_wg1-laynt+1));
    
    
end; % for laynt

% AcaGEx_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% AcaGEy_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% AcaGEz_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% 
% AcaGHx_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% AcaGHy_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% AcaGHz_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);


AcaGEx_xz1=aGEx_xz1(:,:,1);
AcaGEy_xz1=aGEy_xz1(:,:,1);
AcaGEz_xz1=aGEz_xz1(:,:,1);

AcaGHx_xz1=aGHx_xz1(:,:,1);
AcaGHy_xz1=aGHy_xz1(:,:,1);
AcaGHz_xz1=aGHz_xz1(:,:,1);

for laynt=2:SNlay_wg1
    AcaGEx_xz1=[AcaGEx_xz1 aGEx_xz1(:,:,laynt)];
    AcaGEy_xz1=[AcaGEy_xz1 aGEy_xz1(:,:,laynt)];
    AcaGEz_xz1=[AcaGEz_xz1 aGEz_xz1(:,:,laynt)];
    
    AcaGHx_xz1=[AcaGHx_xz1 aGHx_xz1(:,:,laynt)];
    AcaGHy_xz1=[AcaGHy_xz1 aGHy_xz1(:,:,laynt)];
    AcaGHz_xz1=[AcaGHz_xz1 aGHz_xz1(:,:,laynt)];
end;





% AcaGEx_xz2reshape(AcaGEx_xz,size(AcaGEx_xz,1),size(AcaGEx_xz,2)*SNlay_wg1);
% 
% 
% figure(1); imagesc(real([AcaGEx_xz1 ])); colorbar;  title('Ex-xz');
% figure(2); imagesc(real([AcaGEy_xz1 ])); colorbar;  title('Ey-xz');
% figure(3); imagesc(real([AcaGEz_xz1 ])); colorbar;  title('Ez-xz'); 
% 
% figure(4); imagesc(real([AcaGHx_xz1 ])); colorbar;  title('Hx-xz'); 
% figure(5); imagesc(real([AcaGHy_xz1 ])); colorbar;  title('Hy-xz');
% figure(6); imagesc(real([AcaGHz_xz1 ])); colorbar;  title('Hz-xz');
% 
figure(11);
subplot(3,2,1); imagesc((real(AcaGEx_xz1))); title('aPG-Ex-xz'); colorbar;
subplot(3,2,2); imagesc((real(AcaGEy_xz1))); title('aPG-Ey-xz'); colorbar;
subplot(3,2,3); imagesc((real(AcaGEz_xz1))); title('aPG-Ez-xz'); colorbar;
subplot(3,2,4); imagesc((real(AcaGHx_xz1))); title('aPG-Hx-xz'); colorbar;
subplot(3,2,5); imagesc((real(AcaGHy_xz1))); title('aPG-Hy-xz'); colorbar;
subplot(3,2,6); imagesc((real(AcaGHz_xz1))); title('aPG-Hz-xz'); colorbar;



%% field visualization - region III

SNlay_wg3=1;

[xlen zlen]=size(aNG_Ey_xz_wg(:,:,1));  
aGEx_xz3=zeros(xlen,zlen,SNlay_wg3);
aGEy_xz3=zeros(xlen,zlen,SNlay_wg3);
aGEz_xz3=zeros(xlen,zlen,SNlay_wg3);
aGHx_xz3=zeros(xlen,zlen,SNlay_wg3);
aGHy_xz3=zeros(xlen,zlen,SNlay_wg3);
aGHz_xz3=zeros(xlen,zlen,SNlay_wg3);

for laynt=1:SNlay_wg3
    for pbm_ind=1:pcnt
        
        aGEy_xz3(:,:,laynt)=aGEy_xz3(:,:,laynt)+CT3(pbm_ind)*aPG_Ey_xz_wg(:,:,pbm_ind)*exp(ap_evalue_wg(pbm_ind)*aTz_wg*(laynt-1)); %Sx
        aGEx_xz3(:,:,laynt)=aGEx_xz3(:,:,laynt)+CT3(pbm_ind)*aPG_Ex_xz_wg(:,:,pbm_ind)*exp(ap_evalue_wg(pbm_ind)*aTz_wg*(laynt-1)); %Sx
        aGEz_xz3(:,:,laynt)=aGEz_xz3(:,:,laynt)+CT3(pbm_ind)*aPG_Ez_xz_wg(:,:,pbm_ind)*exp(ap_evalue_wg(pbm_ind)*aTz_wg*(laynt-1));
        
        aGHy_xz3(:,:,laynt)=aGHy_xz3(:,:,laynt)+CT3(pbm_ind)*aPG_Hy_xz_wg(:,:,pbm_ind)*exp(ap_evalue_wg(pbm_ind)*aTz_wg*(laynt-1)); %Sx
        aGHx_xz3(:,:,laynt)=aGHx_xz3(:,:,laynt)+CT3(pbm_ind)*aPG_Hx_xz_wg(:,:,pbm_ind)*exp(ap_evalue_wg(pbm_ind)*aTz_wg*(laynt-1)); %Sx
        aGHz_xz3(:,:,laynt)=aGHz_xz3(:,:,laynt)+CT3(pbm_ind)*aPG_Hz_xz_wg(:,:,pbm_ind)*exp(ap_evalue_wg(pbm_ind)*aTz_wg*(laynt-1));
        
    end;
    
    
end; % for laynt

% AcaGEx_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% AcaGEy_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% AcaGEz_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% 
% AcaGHx_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% AcaGHy_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% AcaGHz_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);


AcaGEx_xz3=aGEx_xz3(:,:,1);
AcaGEy_xz3=aGEy_xz3(:,:,1);
AcaGEz_xz3=aGEz_xz3(:,:,1);

AcaGHx_xz3=aGHx_xz3(:,:,1);
AcaGHy_xz3=aGHy_xz3(:,:,1);
AcaGHz_xz3=aGHz_xz3(:,:,1);

for laynt=2:SNlay_wg3
    AcaGEx_xz3=[AcaGEx_xz3 aGEx_xz3(:,:,laynt)];
    AcaGEy_xz3=[AcaGEy_xz3 aGEy_xz3(:,:,laynt)];
    AcaGEz_xz3=[AcaGEz_xz3 aGEz_xz3(:,:,laynt)];
    
    AcaGHx_xz3=[AcaGHx_xz3 aGHx_xz3(:,:,laynt)];
    AcaGHy_xz3=[AcaGHy_xz3 aGHy_xz3(:,:,laynt)];
    AcaGHz_xz3=[AcaGHz_xz3 aGHz_xz3(:,:,laynt)];
end;





% AcaGEx_xz2=reshape(AcaGEx_xz,size(AcaGEx_xz,1),size(AcaGEx_xz,2)*SNlay_wg1);
% 
% 
% figure(1); imagesc(real([AcaGEx_xz1 ])); colorbar;  title('Ex-xz');
% figure(2); imagesc(real([AcaGEy_xz1 ])); colorbar;  title('Ey-xz');
% figure(3); imagesc(real([AcaGEz_xz1 ])); colorbar;  title('Ez-xz'); 
% 
% figure(4); imagesc(real([AcaGHx_xz1 ])); colorbar;  title('Hx-xz'); 
% figure(5); imagesc(real([AcaGHy_xz1 ])); colorbar;  title('Hy-xz');
% figure(6); imagesc(real([AcaGHz_xz1 ])); colorbar;  title('Hz-xz');
% 
figure(13);
subplot(3,2,1); imagesc((real(AcaGEx_xz3))); title('aPG-Ex-xz'); colorbar;
subplot(3,2,2); imagesc((real(AcaGEy_xz3))); title('aPG-Ey-xz'); colorbar;
subplot(3,2,3); imagesc((real(AcaGEz_xz3))); title('aPG-Ez-xz'); colorbar;
subplot(3,2,4); imagesc((real(AcaGHx_xz3))); title('aPG-Hx-xz'); colorbar;
subplot(3,2,5); imagesc((real(AcaGHy_xz3))); title('aPG-Hy-xz'); colorbar;
subplot(3,2,6); imagesc((real(AcaGHz_xz3))); title('aPG-Hz-xz'); colorbar;





%% field visualization - region II

SNlay_resonator=1;

[xlen zlen]=size(aNG_Ey_xz_resonator(:,:,1));  
aGEx_xz2=zeros(xlen,zlen,SNlay_resonator);
aGEy_xz2=zeros(xlen,zlen,SNlay_resonator);
aGEz_xz2=zeros(xlen,zlen,SNlay_resonator);
aGHx_xz2=zeros(xlen,zlen,SNlay_resonator);
aGHy_xz2=zeros(xlen,zlen,SNlay_resonator);
aGHz_xz2=zeros(xlen,zlen,SNlay_resonator);

for laynt=1:SNlay_resonator
    for pbm_ind=1:pcnt
        aGEy_xz2(:,:,laynt)=aGEy_xz2(:,:,laynt)+C_p(pbm_ind,laynt)*aPG_Ey_xz_resonator(:,:,pbm_ind); %Sx
        aGEx_xz2(:,:,laynt)=aGEx_xz2(:,:,laynt)+C_p(pbm_ind,laynt)*aPG_Ex_xz_resonator(:,:,pbm_ind); %Sx
        aGEz_xz2(:,:,laynt)=aGEz_xz2(:,:,laynt)+C_p(pbm_ind,laynt)*aPG_Ez_xz_resonator(:,:,pbm_ind);
        
        aGHy_xz2(:,:,laynt)=aGHy_xz2(:,:,laynt)+C_p(pbm_ind,laynt)*aPG_Hy_xz_resonator(:,:,pbm_ind); %Sx
        aGHx_xz2(:,:,laynt)=aGHx_xz2(:,:,laynt)+C_p(pbm_ind,laynt)*aPG_Hx_xz_resonator(:,:,pbm_ind); %Sx
        aGHz_xz2(:,:,laynt)=aGHz_xz2(:,:,laynt)+C_p(pbm_ind,laynt)*aPG_Hz_xz_resonator(:,:,pbm_ind);
    end;
    for pbm_ind=1:mcnt
        aGEy_xz2(:,:,laynt)=aGEy_xz2(:,:,laynt)+C_n(pbm_ind,laynt)*aNG_Ey_xz_resonator(:,:,pbm_ind); %Sx
        aGEx_xz2(:,:,laynt)=aGEx_xz2(:,:,laynt)+C_n(pbm_ind,laynt)*aNG_Ex_xz_resonator(:,:,pbm_ind); %Sx
        aGEz_xz2(:,:,laynt)=aGEz_xz2(:,:,laynt)+C_n(pbm_ind,laynt)*aNG_Ez_xz_resonator(:,:,pbm_ind);
        
        aGHy_xz2(:,:,laynt)=aGHy_xz2(:,:,laynt)+C_n(pbm_ind,laynt)*aNG_Hy_xz_resonator(:,:,pbm_ind); %Sx
        aGHx_xz2(:,:,laynt)=aGHx_xz2(:,:,laynt)+C_n(pbm_ind,laynt)*aNG_Hx_xz_resonator(:,:,pbm_ind); %Sx
        aGHz_xz2(:,:,laynt)=aGHz_xz2(:,:,laynt)+C_n(pbm_ind,laynt)*aNG_Hz_xz_resonator(:,:,pbm_ind);
    end;
    
    
end; % for laynt

% AcaGEx_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% AcaGEy_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% AcaGEz_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% 
% AcaGHx_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% AcaGHy_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);
% AcaGHz_xz1=zeros(size(aGEx_xz1,1),size(aGEx_xz1,2)*SNlay_wg1);


AcaGEx_xz2=aGEx_xz2(:,:,1);
AcaGEy_xz2=aGEy_xz2(:,:,1);
AcaGEz_xz2=aGEz_xz2(:,:,1);

AcaGHx_xz2=aGHx_xz2(:,:,1);
AcaGHy_xz2=aGHy_xz2(:,:,1);
AcaGHz_xz2=aGHz_xz2(:,:,1);

for laynt=2:SNlay_resonator
    AcaGEx_xz2=[AcaGEx_xz2 aGEx_xz2(:,:,laynt)];
    AcaGEy_xz2=[AcaGEy_xz2 aGEy_xz2(:,:,laynt)];
    AcaGEz_xz2=[AcaGEz_xz2 aGEz_xz2(:,:,laynt)];
    
    AcaGHx_xz2=[AcaGHx_xz2 aGHx_xz2(:,:,laynt)];
    AcaGHy_xz2=[AcaGHy_xz2 aGHy_xz2(:,:,laynt)];
    AcaGHz_xz2=[AcaGHz_xz2 aGHz_xz2(:,:,laynt)];
end;





% AcaGEx_xz2=reshape(AcaGEx_xz,size(AcaGEx_xz,1),size(AcaGEx_xz,2)*SNlay_wg1);
% 
% 
% figure(1); imagesc(real([AcaGEx_xz1 ])); colorbar;  title('Ex-xz');
% figure(2); imagesc(real([AcaGEy_xz1 ])); colorbar;  title('Ey-xz');
% figure(3); imagesc(real([AcaGEz_xz1 ])); colorbar;  title('Ez-xz'); 
% 
% figure(4); imagesc(real([AcaGHx_xz1 ])); colorbar;  title('Hx-xz'); 
% figure(5); imagesc(real([AcaGHy_xz1 ])); colorbar;  title('Hy-xz');
% figure(6); imagesc(real([AcaGHz_xz1 ])); colorbar;  title('Hz-xz');
% 
figure(12);
subplot(3,2,1); imagesc((real(AcaGEx_xz2))); title('aPG-Ex-xz'); colorbar;
subplot(3,2,2); imagesc((real(AcaGEy_xz2))); title('aPG-Ey-xz'); colorbar;
subplot(3,2,3); imagesc((real(AcaGEz_xz2))); title('aPG-Ez-xz'); colorbar;
subplot(3,2,4); imagesc((real(AcaGHx_xz2))); title('aPG-Hx-xz'); colorbar;
subplot(3,2,5); imagesc((real(AcaGHy_xz2))); title('aPG-Hy-xz'); colorbar;
subplot(3,2,6); imagesc((real(AcaGHz_xz2))); title('aPG-Hz-xz'); colorbar;




%% hap

AcaGEx_xz_total=[AcaGEx_xz1 AcaGEx_xz2 AcaGEx_xz3];
AcaGEy_xz_total=[AcaGEy_xz1 AcaGEy_xz2 AcaGEy_xz3];
AcaGEz_xz_total=[AcaGEz_xz1 AcaGEz_xz2 AcaGEz_xz3];
AcaGHx_xz_total=[AcaGHx_xz1 AcaGHx_xz2 AcaGHx_xz3];
AcaGHy_xz_total=[AcaGHy_xz1 AcaGHy_xz2 AcaGHy_xz3];
AcaGHz_xz_total=[AcaGHz_xz1 AcaGHz_xz2 AcaGHz_xz3];


figure(1);
subplot(3,2,1); imagesc((real(AcaGEx_xz_total))); title('aPG-Ex-xz'); colorbar;
subplot(3,2,2); imagesc((real(AcaGEy_xz_total))); title('aPG-Ey-xz'); colorbar;
subplot(3,2,3); imagesc((real(AcaGEz_xz_total))); title('aPG-Ez-xz'); colorbar;
subplot(3,2,4); imagesc((real(AcaGHx_xz_total))); title('aPG-Hx-xz'); colorbar;
subplot(3,2,5); imagesc((real(AcaGHy_xz_total))); title('aPG-Hy-xz'); colorbar;
subplot(3,2,6); imagesc((real(AcaGHz_xz_total))); title('aPG-Hz-xz'); colorbar;


