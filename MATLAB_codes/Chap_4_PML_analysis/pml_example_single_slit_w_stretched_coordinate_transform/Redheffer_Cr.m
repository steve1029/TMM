function [Ca Cb]=Redheffer_Cr(R1a,T1a,R1b,T1b,R2a,T2a,R2b,T2b,inv2a1b,inv1b2a,Ca_prev,Cb_prev)

Ca=          Ca_prev*inv1b2a*T1a;       % eq. 3.2.20c on page 36
Cb=Cb_prev + Ca_prev*inv1b2a*R1b*T2b;   % eq. 3.2.20d on page 36

