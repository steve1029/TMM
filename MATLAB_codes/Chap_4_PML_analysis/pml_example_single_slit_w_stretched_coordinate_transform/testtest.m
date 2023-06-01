clear all
close all
clc
%%
L=50;
Nlay=2;

%%
RRa=rand(2*L,2*L);
TTa=rand(2*L,2*L);
RRb=rand(2*L,2*L);
TTb=rand(2*L,2*L);
Lfree_Rf=rand(2*L,2*L);
Lfree_Tf=rand(2*L,2*L);
Lfree_Rb=rand(2*L,2*L);
Lfree_Tb=rand(2*L,2*L);
Rwg_Tf2=rand(2*L,2*L);
Rwg_Rb2=rand(2*L,2*L);
Rwg_Tb2=rand(2*L,2*L);
Rwg_Rf2=rand(2*L,2*L);
Ca=rand(2*L,2*L,Nlay);
Cb=rand(2*L,2*L,Nlay);

save rand_result;

%%
tt=clock;
Bdr_Smat_case2;
t_v1=etime(clock,tt);
m1=RRb;

%%
load rand_result;
tt=clock;
Bdr_Smat_case2_v2;
t_v2=etime(clock,tt);
m2=RRb;

%%
load rand_result;
tt=clock;
Bdr_Smat_case2_v3;
t_v3=etime(clock,tt);
m3=RRb;

%%
disp('t1  t2  t3');
disp([num2str(t_v1) ' ' num2str(t_v2) ' ' num2str(t_v3)]);

disp(['diff12 :  ' num2str(sum(sum(abs(m1-m2))))]);
disp(['diff23 :  ' num2str(sum(sum(abs(m2-m3))))]);
disp(['diff13 :  ' num2str(sum(sum(abs(m1-m3))))]);




