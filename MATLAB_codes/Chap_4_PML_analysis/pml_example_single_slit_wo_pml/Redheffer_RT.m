function [Ra Ta Rb Tb]=Redheffer_RT(R1a,T1a,R1b,T1b,R2a,T2a,R2b,T2b,inv2a1b,inv1b2a)

Ra=R1a + T1b*inv2a1b*R2a*T1a;  % eq. 3.2.17a on page 34
Ta=      T2a*inv1b2a*T1a;      % eq. 3.2.17b on page 34
Rb=R2b + T2a*inv1b2a*R1b*T2b;  % eq. 3.2.17c on page 34
Tb=      T1b*inv2a1b*T2b;      % eq. 3.2.17d on page 34
