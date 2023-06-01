clear all
% close all
clc

%%
% 제1의 자유도는 입사각도.
% 그 외에는 두께. 파장 정보는 두께에 대한 상대적 비율로 들어가게 된다.
% TE 나 TM 인 경우도 있겠구나.
% 우선 계단 모양으로 급격하게 꺾이는 경우를 보자.

% param_txt       = ['mu'];
% param_set       = linspace(-1.001,-0.999,1001);param_len=length(param_set);

param_txt       = 'angle';
param_set       = linspace(0,89,90);param_len=length(param_set);

% param_txt       = ['\tau'];
% param_set       = linspace(0.8,1.2,1001);param_len=length(param_set);


%% 계산
answer_t        = zeros(1,param_len);
answer_r        = zeros(1,param_len);

for param_idx=1:param_len
    param       = param_set(param_idx);
    
    lambda      = 532e-9;k0=2*pi/lambda;w=3e8*k0;
    inc_ang     = param*pi/180;
    mode        = 1; % TM
    ep          = [1 1+i 1].^2;
    mu          = [1 1 1];
    d           = [0 1 0]*4*lambda;N=length(ep);

%     simple_tmm_callee;
    simple_etmm_callee_inner_coeff;
    
    answer_t(param_idx) = trans*conj(trans);
    answer_r(param_idx) = refle*conj(refle);
end 


%% 결과 표시
% close all
sum_d=zeros(1,length(d));
for k=1:length(d)
    sum_d(k)=sum(d(1:k))/lambda;
end
    
figure(11);set(gca,'fontsize',16);set(gca,'fontname','times new roman');hold on;grid on;box on;
plot(sum_d(2:end-1)*linspace(0,5,101),real(sqrt(ep(2:end-1)))*ones(1,101),'-r','linewidth',3);
plot(sum_d(2:end-1)*linspace(0,5,101),imag(sqrt(ep(2:end-1)))*ones(1,101),':b','linewidth',3);
axis([0 5.2 0 2]);set(gca,'fontname','times new roman');
legend('Re(n_A_B_L)','Im(n_A_B_L)');

figure(12);set(gca,'fontsize',16);set(gca,'fontname','times new roman');hold on;grid on;box on;
plot(param_set,answer_r,'-r','linewidth',3);
plot(param_set,answer_t,':b','linewidth',3);
axis([param_set(1) param_set(end) 0 1]);set(gca,'fontname','times new roman');
% xlabel('Incident angle (deg.)');set(gca,'fontname','times new roman');
legend('Reflection','Transmission');

figure(13);set(gca,'fontsize',16);set(gca,'fontname','times new roman');hold on;grid on;box on;
plot(param_set,10*log10(answer_r),'-r','linewidth',3);
plot(param_set,10*log10(answer_t),':b','linewidth',3);
axis([param_set(1) param_set(end) -300 0]);set(gca,'fontname','times new roman');
% xlabel('Incident angle (deg.)');set(gca,'fontname','times new roman');
legend('Reflection','Transmission');
