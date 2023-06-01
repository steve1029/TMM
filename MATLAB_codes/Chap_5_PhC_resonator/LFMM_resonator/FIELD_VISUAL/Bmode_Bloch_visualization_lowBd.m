% Bmode field visualization in x-z plane
% Lower grating
% freespace
%    |
%  grating

D=zeros(2*L,1);
inc_mode=22;
D(inc_mode)=1;


C_n=BLR33*D;        % negative coupling coefficient
T4_temp=BLT34*D;    % positive coupling coefficient


T4_z=zeros(L,1);
T4_y=zeros(L,1);
T4_y=T4_temp(1:L);
T4_z=T4_temp(L+1:2*L);


T4x_cof=zeros(NBz,NBy);
T4y_cof=zeros(NBz,NBy);
T4z_cof=zeros(NBz,NBy);


for k=1:NBz
   for l=1:NBy
      
      T4y_cof(k,l)=T4_y((k-1)*NBy+l);
      T4z_cof(k,l)=T4_z((k-1)*NBy+l );
      T4x_cof(k,l)=-( Bkz_vc(k)*T4z_cof(k,l)+Bky_vc(l)*T4y_cof(k,l) )/Bkx_vc(k,l); 

   end;
end;

% region 4
Ey_xz4=zeros(length(xx4),length(zz));
Ex_xz4=zeros(length(xx4),length(zz));
Ez_xz4=zeros(length(xx4),length(zz));

[z4 x4]=meshgrid(zz,xx4); % region 4

 for k=1:NBz
         for l=1:NBy
           Ey_xz4=Ey_xz4+T4y_cof(k,l)*exp(j*( Bkz_vc(k)*z4+Bky_vc(l)*y+Bkx_vc(k,l)*x4 ));   
           Ex_xz4=Ex_xz4+T4x_cof(k,l)*exp(j*( Bkz_vc(k)*z4+Bky_vc(l)*y+Bkx_vc(k,l)*x4 ));
           Ez_xz4=Ez_xz4+T4z_cof(k,l)*exp(j*( Bkz_vc(k)*z4+Bky_vc(l)*y+Bkx_vc(k,l)*x4 ));
           end;
        end;


% Grating region
      
bGEx_xz=zeros(3*length(xx),length(zz));
bGEy_xz=zeros(3*length(xx),length(zz));
bGEz_xz=zeros(3*length(xx),length(zz));


   for pbm_ind=1:mcnt
                        
            bGEy_xz=bGEy_xz+C_n(pbm_ind)*[bNG_Ey_xz(:,:,pbm_ind)*exp(-bm_evalue(pbm_ind)*2*Tx);  bNG_Ey_xz(:,:,pbm_ind)*exp(-bm_evalue(pbm_ind)*Tx);   bNG_Ey_xz(:,:,pbm_ind)]; 
            bGEx_xz=bGEx_xz+C_n(pbm_ind)*[bNG_Ex_xz(:,:,pbm_ind)*exp(-bm_evalue(pbm_ind)*2*Tx) ; bNG_Ex_xz(:,:,pbm_ind)*exp(-bm_evalue(pbm_ind)*Tx);   bNG_Ex_xz(:,:,pbm_ind)]; 
         	bGEz_xz=bGEz_xz+C_n(pbm_ind)*[bNG_Ez_xz(:,:,pbm_ind)*exp(-bm_evalue(pbm_ind)*2*Tx);  bNG_Ez_xz(:,:,pbm_ind)*exp(-bm_evalue(pbm_ind)*Tx) ;  bNG_Ez_xz(:,:,pbm_ind)]; 
                  
    end;
    
            bGEy_xz=bGEy_xz+[bPG_Ey_xz(:,:,inc_mode)*exp(-bp_evalue(inc_mode)*2*Tx) ;bPG_Ey_xz(:,:,inc_mode)*exp(-bp_evalue(inc_mode)*Tx) ; bPG_Ey_xz(:,:,inc_mode)]; 
            bGEx_xz=bGEx_xz+[bPG_Ex_xz(:,:,inc_mode)*exp(-bp_evalue(inc_mode)*2*Tx); bPG_Ex_xz(:,:,inc_mode)*exp(-bp_evalue(inc_mode)*Tx) ; bPG_Ex_xz(:,:,inc_mode)]; 
          	bGEz_xz=bGEz_xz+[bPG_Ez_xz(:,:,inc_mode)*exp(-bp_evalue(inc_mode)*2*Tx); bPG_Ez_xz(:,:,inc_mode)*exp(-bp_evalue(inc_mode)*Tx) ; bPG_Ez_xz(:,:,inc_mode)]; 
%        
   
figure(1); mesh(100*real([bGEx_xz; Ex_xz4])); title('Ex-xz'); title('Ex-xz');axis equal; 
figure(2); mesh(100*real([bGEy_xz; Ey_xz4])); title('Ey-xz'); title('Ey-xz');axis equal; 
figure(3); mesh(100*real([bGEz_xz; Ez_xz4])); title('Ez-xz'); title('Ez-xz');axis equal; 




