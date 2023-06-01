% Bmode field visualization in x-z plane
% Lower grating
% grating
%    |
% freespace

D=zeros(2*L,1);
d11=zeros(L,1);
d12=zeros(L,1);
centx=nx*NBy+ny+1;
centy=nx*NBy+ny+1;
d11(centx)=Uy;  % Kx of incident beam
d12(centy)=Ux;
D(1:L)=d11;
D(L+1:2*L)=d12;
% D=zeros(2*L,1);
% inc_mode=22;
% D(inc_mode)=1;

T3_temp=BRR33*D;
C_p=BRT34*D;

T3_z=zeros(L,1);
T3_y=zeros(L,1);
T3_y=T3_temp(1:L);
T3_z=T3_temp(L+1:2*L);


T3x_cof=zeros(NBz,NBy);
T3y_cof=zeros(NBz,NBy);
T3z_cof=zeros(NBz,NBy);


for k=1:NBz
   for l=1:NBy
      
      T3y_cof(k,l)=T3_y((k-1)*NBy+l);
      T3z_cof(k,l)=T3_z((k-1)*NBy+l );
      T3x_cof(k,l)=( Bkz_vc(k)*T3z_cof(k,l)+Bky_vc(l)*T3y_cof(k,l) )/Bkx_vc(k,l); 

   end;
end;

% region 3
Ey_xz3=zeros(length(xx3),length(zz));
Ex_xz3=zeros(length(xx3),length(zz));
Ez_xz3=zeros(length(xx3),length(zz));

[z3 x3]=meshgrid(zz,xx3); % region 4
y=0;
 for k=1:NBz
         for l=1:NBy
           Ey_xz3=Ey_xz3+T3y_cof(k,l)*exp(j*( Bkz_vc(k)*z3+Bky_vc(l)*y-Bkx_vc(k,l)*x3 ));   
           Ex_xz3=Ex_xz3+T3x_cof(k,l)*exp(j*( Bkz_vc(k)*z3+Bky_vc(l)*y-Bkx_vc(k,l)*x3 ));
           Ez_xz3=Ez_xz3+T3z_cof(k,l)*exp(j*( Bkz_vc(k)*z3+Bky_vc(l)*y-Bkx_vc(k,l)*x3 ));
           end;
        end;
       Ey_xz3=Ey_xz3+Uy*exp(j*( kx0*z3+ky0*y+kz0*x3));
       Ex_xz3=Ex_xz3+Uz*exp(j*( kx0*z3+ky0*y+kz0*x3));
       Ez_xz3=Ez_xz3+Ux*exp(j*( kx0*z3+ky0*y+kz0*x3));

% Grating region
      
bGEx_xz=zeros(3*length(xx),length(zz));
bGEy_xz=zeros(3*length(xx),length(zz));
bGEz_xz=zeros(3*length(xx),length(zz));


   for pbm_ind=1:pcnt
                        
            bGEy_xz=bGEy_xz+C_p(pbm_ind)*[bPG_Ey_xz(:,:,pbm_ind); bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*Tx); bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*Tx)];  
            bGEx_xz=bGEx_xz+C_p(pbm_ind)*[bPG_Ex_xz(:,:,pbm_ind); bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*Tx); bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*Tx)]; 
         	bGEz_xz=bGEz_xz+C_p(pbm_ind)*[bPG_Ez_xz(:,:,pbm_ind); bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*Tx); bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*Tx)]; 
                  
    end;
    

figure(1); mesh(100*abs([Ex_xz3; bGEx_xz])); title('Ex-xz');axis equal; 
figure(2); mesh(100*abs([Ey_xz3; bGEy_xz]));  title('Ey-xz');axis equal; 
figure(3); mesh(100*abs([Ez_xz3; bGEz_xz])); title('Ez-xz');axis equal; 


