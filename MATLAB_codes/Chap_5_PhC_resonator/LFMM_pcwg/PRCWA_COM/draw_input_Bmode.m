
%% input mode visualization

for inc_mode=1:pcnt
%inc_mode=1;
pbm_ind=inc_mode;

% E-field

bGEy_xz=[bPG_Ey_xz(:,:,pbm_ind);  
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx) ; 
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*bTx);
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*3*bTx) ;
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*4*bTx);
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*5*bTx);
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*6*bTx);  
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*7*bTx)]; 

bGEx_xz=[bPG_Ex_xz(:,:,pbm_ind);  
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx);
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*bTx);
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*3*bTx);
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*4*bTx);
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*5*bTx);
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*6*bTx);
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*7*bTx)]; 

bGEz_xz=[bPG_Ez_xz(:,:,pbm_ind);
    bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx);
    bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*bTx);
      bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*3*bTx);
    bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*4*bTx);
      bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*5*bTx);
    bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*6*bTx);
     bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*7*bTx);
    ]; 


% H-field

bGHy_xz=[bPG_Hy_xz(:,:,pbm_ind); 
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx); 
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*bTx);
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*3*bTx); 
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*4*bTx);
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*5*bTx);
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*6*bTx); 
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*7*bTx)]; 

bGHx_xz=[bPG_Hx_xz(:,:,pbm_ind);  
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx);
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*bTx);
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*3*bTx);
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*4*bTx);
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*5*bTx);
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*6*bTx);
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*7*bTx)]; 

bGHz_xz=[bPG_Hz_xz(:,:,pbm_ind);
    bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx);
    bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*bTx);
      bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*3*bTx);
    bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*4*bTx);
      bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*5*bTx);
    bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*6*bTx);
     bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*7*bTx)]; 
                  
figure(11); imagesc(abs(real(bGEx_xz))); colorbar;
figure(12); imagesc(abs(real(bGEy_xz))); colorbar;
figure(13); imagesc(abs(real(bGEz_xz))); colorbar;

figure(4); imagesc(abs(real(bGHx_xz))); colorbar;
figure(5); imagesc(abs(real(bGHy_xz))); colorbar;
figure(6); imagesc(abs(real(bGHz_xz))); colorbar;

pause;
end;
% 
% [lx lz]=size(aGEx_xz);
% 
% axx=xx;
% x_inc=axx(2)-axx(1);
% azz=[0:z_inc:z_inc*(lz-1)];
% 
% 
% figure(1); imagesc(azz/a,axx/a,100*real(aGEy_xz));
% axis([azz(1)/a azz(lz)/a axx(1)/a axx(lx)/a ]);
% hold on;
% 
% for m=1:8
%     
%     for n=-4:4
%     
%     rectangle('Position',[m-r/a-0.5,n-r/a,2*r/a,2*r/a],'Curvature',[1,1],'LineWidth',2,'LineStyle','--','EdgeColor','w');
%         
%     end;
%     
% end;
% 
% 
for inc_mode=1:pcnt
%inc_mode=1;
pbm_ind=inc_mode;

% E-field

bGEy_xz=[bPG_Ey_xz(:,:,pbm_ind);  
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx) ; 
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*bTx);
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*3*bTx) ;
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*4*bTx);
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*5*bTx);
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*6*bTx);  
    bPG_Ey_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*7*bTx)]; 

bGEx_xz=[bPG_Ex_xz(:,:,pbm_ind);  
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx);
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*bTx);
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*3*bTx);
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*4*bTx);
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*5*bTx);
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*6*bTx);
    bPG_Ex_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*7*bTx)]; 

bGEz_xz=[bPG_Ez_xz(:,:,pbm_ind);
    bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx);
    bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*bTx);
      bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*3*bTx);
    bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*4*bTx);
      bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*5*bTx);
    bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*6*bTx);
     bPG_Ez_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*7*bTx);
    ]; 


% H-field

bGHy_xz=[bPG_Hy_xz(:,:,pbm_ind); 
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx); 
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*bTx);
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*3*bTx); 
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*4*bTx);
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*5*bTx);
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*6*bTx); 
    bPG_Hy_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*7*bTx)]; 

bGHx_xz=[bPG_Hx_xz(:,:,pbm_ind);  
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx);
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*bTx);
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*3*bTx);
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*4*bTx);
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*5*bTx);
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*6*bTx);
    bPG_Hx_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*7*bTx)]; 

bGHz_xz=[bPG_Hz_xz(:,:,pbm_ind);
    bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*bTx);
    bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*2*bTx);
      bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*3*bTx);
    bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*4*bTx);
      bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*5*bTx);
    bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*6*bTx);
     bPG_Hz_xz(:,:,pbm_ind)*exp(bp_evalue(pbm_ind)*7*bTx)]; 
                  
figure(11); imagesc(abs(real(bGEx_xz))); colorbar;
figure(12); imagesc(abs(real(bGEy_xz))); colorbar;
figure(13); imagesc(abs(real(bGEz_xz))); colorbar;

figure(4); imagesc(abs(real(bGHx_xz))); colorbar;
figure(5); imagesc(abs(real(bGHy_xz))); colorbar;
figure(6); imagesc(abs(real(bGHz_xz))); colorbar;

pause;
end;