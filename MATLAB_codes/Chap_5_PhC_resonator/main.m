clear all
close all
clc

%% run
cd LFMM_pcwg;LFMM_pcwg;cd ..;
cd LFMM_resonator;LFMM_resonator;cd ..;


%% pcwg
cd LFMM_pcwg;
load 2portDATA_pcwg;
cd ..

% left waveguide
Lwg_Tf=ALT12;   clear ALT12;
Lwg_Rb=ALR11;   clear ALR11;
Lwg_Tb=ALT21;   clear ALT21;
Lwg_Rf=ALR22;   clear ALR22;

% right waveguide
Rwg_Tf=ART12;  clear ART12;
Rwg_Rb=ARR11;  clear ARR11;
Rwg_Tb=ART21;  clear ART21;
Rwg_Rf=ARR22;  clear ARR22;

% eigenmode field distribution
aPG_Ex_xz_wg=aPG_Ex_xz;  clear aPG_Ex_xz;
aPG_Ey_xz_wg=aPG_Ey_xz;  clear aPG_Ey_xz;
aPG_Ez_xz_wg=aPG_Ez_xz;  clear aPG_Ez_xz;
aPG_Hx_xz_wg=aPG_Hx_xz;  clear aPG_Hx_xz;
aPG_Hy_xz_wg=aPG_Hy_xz;  clear aPG_Hy_xz;
aPG_Hz_xz_wg=aPG_Hz_xz;  clear aPG_Hz_xz;
aNG_Ex_xz_wg=aNG_Ex_xz;  clear aNG_Ex_xz;
aNG_Ey_xz_wg=aNG_Ey_xz;  clear aNG_Ey_xz;
aNG_Ez_xz_wg=aNG_Ez_xz;  clear aNG_Ez_xz;
aNG_Hx_xz_wg=aNG_Hx_xz;  clear aNG_Hx_xz;
aNG_Hy_xz_wg=aNG_Hy_xz;  clear aNG_Hy_xz;
aNG_Hz_xz_wg=aNG_Hz_xz;  clear aNG_Hz_xz;

% eigenvalue
ap_evalue_wg=ap_evalue; clear ap_evalue;
am_evalue_wg=am_evalue; clear am_evalue;

aTz_wg=aTz;



%% resonator

cd LFMM_resonator;
load 2portDATA_resonator;
cd ..

% layer S matrix
RRa_resonator=AR11;  clear AR11;
TTa_resonator=AT12;  clear AT12;
RRb_resonator=AR22;  clear AR22;
TTb_resonator=AT21;  clear AT21;

% eigenmode field distribution
aPG_Ex_xz_resonator=aPG_Ex_xz;  clear aPG_Ex_xz;
aPG_Ey_xz_resonator=aPG_Ey_xz;  clear aPG_Ey_xz;
aPG_Ez_xz_resonator=aPG_Ez_xz;  clear aPG_Ez_xz;
aPG_Hx_xz_resonator=aPG_Hx_xz;  clear aPG_Hx_xz;
aPG_Hy_xz_resonator=aPG_Hy_xz;  clear aPG_Hy_xz;
aPG_Hz_xz_resonator=aPG_Hz_xz;  clear aPG_Hz_xz;
aNG_Ex_xz_resonator=aNG_Ex_xz;  clear aNG_Ex_xz;
aNG_Ey_xz_resonator=aNG_Ey_xz;  clear aNG_Ey_xz;
aNG_Ez_xz_resonator=aNG_Ez_xz;  clear aNG_Ez_xz;
aNG_Hx_xz_resonator=aNG_Hx_xz;  clear aNG_Hx_xz;
aNG_Hy_xz_resonator=aNG_Hy_xz;  clear aNG_Hy_xz;
aNG_Hz_xz_resonator=aNG_Hz_xz;  clear aNG_Hz_xz;

% eigenvalue
ap_evalue_resonator=ap_evalue; clear ap_evalue;
am_evalue_resonator=am_evalue; clear am_evalue;

ACa_resonator=ACa;
ACb_resonator=ACb;

aTz_resonator=aTz;

%% Final boundary matching
%% .    input + body
%
I=eye(2*L,2*L);
T_temp1a=Lwg_Tf;       %outTb
R_temp1a=Lwg_Rb;       %outRf
T_temp1b=Lwg_Tb;       %outTf
R_temp1b=Lwg_Rf;       %outRb

%%% Important
     Ca=ACa_resonator;
     Cb=ACb_resonator;
%%%
%
T_temp2a=TTa_resonator;
R_temp2a=RRa_resonator;
T_temp2b=TTb_resonator;
R_temp2b=RRb_resonator;

RRa=(R_temp1a+T_temp1b*inv(I-R_temp2a*R_temp1b)*R_temp2a*T_temp1a);
TTa=T_temp2a*inv(I-R_temp1b*R_temp2a)*T_temp1a;
%
RRb=(R_temp2b+T_temp2a*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b);
TTb=T_temp1b*inv(I-R_temp2a*R_temp1b)*T_temp2b;

clear tCa;
clear tCb;

for k=1:1
        tCa(:,:,k)=Ca(:,:,k)*inv(I-R_temp1b*R_temp2a)*T_temp1a;
        tCb(:,:,k)=Cb(:,:,k)+Ca(:,:,k)*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b;
end;
    Ca=tCa;
    Cb=tCb;

%% .    (input + body) + output
%
I=eye(2*L,2*L);
T_temp1a=TTa;       %outTb
R_temp1a=RRa;       %outRf
T_temp1b=TTb;       %outTf
R_temp1b=RRb;       %outRb

T_temp2a=Rwg_Tf;
R_temp2a=Rwg_Rb;
T_temp2b=Rwg_Tb;
R_temp2b=Rwg_Rf;

RRa=(R_temp1a+T_temp1b*inv(I-R_temp2a*R_temp1b)*R_temp2a*T_temp1a);
TTa=T_temp2a*inv(I-R_temp1b*R_temp2a)*T_temp1a;
%
RRb=(R_temp2b+T_temp2a*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b);
TTb=T_temp1b*inv(I-R_temp2a*R_temp1b)*T_temp2b;

for k=1:1
    tCa(:,:,k)=Ca(:,:,k)+Cb(:,:,k)*inv(I-R_temp2a*R_temp1b)*R_temp2a*T_temp1a;
    tCb(:,:,k)=Cb(:,:,k)*inv(I-R_temp2a*R_temp1b)*T_temp2b;
end;
    Ca=tCa;
    Cb=tCb;

%% Field visualization

close all;
LFMM_Amode_field_visual_xz_case4_Lwg_Rwg_leftright;
%LFMM_Amode_field_visual_xz_case4_Lwg_Rwg_rightleft;
