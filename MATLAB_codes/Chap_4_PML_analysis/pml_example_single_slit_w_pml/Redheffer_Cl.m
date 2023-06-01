function [Ca Cb]=Redheffer_Cl(R1a,T1a,R1b,T1b,R2a,T2a,R2b,T2b,inv2a1b,inv1b2a,Ca_prev,Cb_prev)

Ca=Ca_prev + Cb_prev*inv2a1b*R2a*T1a;   % eq. 3.2.20a on page 36
Cb=          Cb_prev*inv2a1b*T2b;       % eq. 3.2.20b on page 36

