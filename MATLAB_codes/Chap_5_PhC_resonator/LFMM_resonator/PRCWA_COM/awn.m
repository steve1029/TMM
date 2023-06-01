function [ey ex]=awn(zv,k,l,pbm_ind,EFy,EFx,me_value)


global aTx; global aTy; global aTz; global nx; global ny; global nz; 
global NBx; global NBy; global NBz; global num_hx; global num_hy; 

Tx=aTx;
Ty=aTy;
Tz=aTz;

ey=0;
ex=0;

for m=1:NBz
  ey=ey+EFy(k,l,m,pbm_ind)*exp(j*(2*pi/Tz*(m-nz-1)*zv))*exp(me_value(pbm_ind)*zv); %Sx
  ex=ex+EFx(k,l,m,pbm_ind)*exp(j*(2*pi/Tz*(m-nz-1)*zv))*exp(me_value(pbm_ind)*zv); %Sx
end;