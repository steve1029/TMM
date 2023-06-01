function [ey ez]=byn(xv,k,l,pbm_ind,EFy,EFz,me_value)


global bTx; global bTy; global bTz; global nx; global ny; global nz; 
global NBx; global NBy; global NBz; global num_hx; global num_hy; 


Tx=bTx;
Ty=bTy;
Tz=bTz;

ey=0;
ez=0;

for m=1:NBz
  ey=ey+EFy(m,l,k,pbm_ind)*exp(j*( (2*pi/Tx*(m-nz-1))*xv))*exp(me_value(pbm_ind)*xv); %Sx
  ez=ez+EFz(m,l,k,pbm_ind)*exp(j*( (2*pi/Tx*(m-nz-1))*xv))*exp(me_value(pbm_ind)*xv); %Sx
end;
            