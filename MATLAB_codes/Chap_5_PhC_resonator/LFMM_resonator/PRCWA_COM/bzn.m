function [hy hz]=bzn(xv,k,l,pbm_ind,HFy,HFz,me_value)


global bTx; global bTy; global bTz; global nx; global ny; global nz; 
global NBx; global NBy; global NBz; global num_hx; global num_hy; 

Tx=bTx;
Ty=bTy;
Tz=bTz;

hy=0;
hz=0;

for m=1:NBz
  hy=hy+HFy(m,l,k,pbm_ind)*exp(j*( (2*pi/Tx*(m-nz-1))*xv ))*exp(me_value(pbm_ind)*xv); %Sx
  hz=hz+HFz(m,l,k,pbm_ind)*exp(j*( (2*pi/Tx*(m-nz-1))*xv ))*exp(me_value(pbm_ind)*xv); %Sx
end;
            