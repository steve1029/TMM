
 
%%%%% Mode power of A_mode 

Iden=zeros(2*L,2*L);
for k=1:2*L
Iden(k,k)=1;   
end;
      
% L+1 layer : region II
KII=zeros(2*L,2*L); % KII matrix
for k=1:NBx
   for l=1:NBy
od_ind1=(k-1)*NBy+l;      
K11(od_ind1,od_ind1)=(kx_vc(k)*ky_vc(l))/kz_vc(k,l);
K12(od_ind1,od_ind1)=(kz_vc(k,l)^2+kx_vc(k)^2)/kz_vc(k,l);
K21(od_ind1,od_ind1)=-(ky_vc(l)^2+kz_vc(k,l)^2)/kz_vc(k,l);
K22(od_ind1,od_ind1)=-( ky_vc(l)*kx_vc(k) )/kz_vc(k,l);

   end;
end;

KII(1:L, 1:L)=K11;
KII(1:L, L+1:2*L)=K12;
KII(L+1:2*L, 1:L)=K21;
KII(L+1:2*L,L+1:2*L )=K22;

Wh=Iden;
Vh=KII/(w0*mu0);

PS=[    Wh                  Wh
        Vh                 -Vh];

%%% power calculation 
bPw=zeros(2*L,1); % z-direction power
bNw=zeros(2*L,1); % z-direction power
% positive_mode
for pbm_ind=1:2*L
    
        tt=PS*bp_evector(:,pbm_ind);
        
        efy=tt(1:L);
        efx=tt(L+1:2*L);
        hfy=tt(2*L+1:3*L);
        hfx=tt(3*L+1:4*L);
    
        bPw(pbm_ind)=0.5*real( sum(efx.*conj(hfy)-efy.*conj(hfx)) );
        
end;
% 
% for pbm_ind=1:2*L
%     pe_vector(:,pbm_ind)=pe_vector(:,pbm_ind)/aPw(pbm_ind)^0.5;
% end;


% negative mode
for pbm_ind=1:2*L
       tt=PS*bm_evector(:,pbm_ind);
        
        efy=tt(1:L);
        efx=tt(L+1:2*L);
        hfy=tt(2*L+1:3*L);
        hfx=tt(3*L+1:4*L);
    
        bNw(pbm_ind)=0.5*real( sum(efx.*conj(hfy)-efy.*conj(hfx)) );
    
end;
% 
% for pbm_ind=1:2*L
%     me_vector(:,pbm_ind)=me_vector(:,pbm_ind)/aNw(pbm_ind)^0.5;
% end;

figure(3); plot(abs(bPw));title('B_positive_mode power');
figure(4); plot(abs(bNw));title('B_negative_mode power');
