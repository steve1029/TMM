function [hy hx]=avp(zv,k,l,pbm_ind,HFy,HFx,pe_value)


global aTx; global aTy; global aTz; global nx; global ny; global nz; 
global NBx; global NBy; global NBz; global num_hx; global num_hy; 

Tx=aTx;
Ty=aTy;
Tz=aTz;

hy=0;
hx=0;

for m=1:NBz
  hy=hy+HFy(k,l,m,pbm_ind)*exp(j*(2*pi/Tz*(m-nz-1)*zv))*exp(pe_value(pbm_ind)*zv);  %Sx
  hx=hx+HFx(k,l,m,pbm_ind)*exp(j*(2*pi/Tz*(m-nz-1)*zv))*exp(pe_value(pbm_ind)*zv);  %Sx
end;
            