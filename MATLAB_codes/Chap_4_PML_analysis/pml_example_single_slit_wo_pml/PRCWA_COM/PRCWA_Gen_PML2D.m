% PRCWA_Gen_PML2D
um=micro;
nm=nano;
%%% PML setting %%%
%-----------------------------------------------------%
pml_width=5*lambda;   % PML thickness
pml_N=30;         % number of staircases along a single direction
%-----------------------------------------------------%

n_surr=n0;
epra=n_surr^2;
% PML basic parameters
pml_num=[1:pml_N];
p_factor=[0.947 1.043 4.552 7.343];
pml_index=n_surr+(p_factor(1)*(pml_num/pml_N)).^p_factor(3)...
    +j*(p_factor(2)*(pml_num/pml_N)).^p_factor(4);

gambda_x=Tx;
pml_epsr=pml_index.^2;
pml_fill_factor=(gambda_x-2*pml_width*fliplr(pml_num)/pml_N)/gambda_x;


Epsr_PML=zeros(num_hx,num_hy);
Apsr_PML=zeros(num_hx,num_hy);

Epsr_nonPML=zeros(num_hx,num_hy);
Apsr_nonPML=zeros(num_hx,num_hy);

f_pml=pml_fill_factor(1);

% non-PML area 
for mm=1:num_hx
    for nn=1:num_hy
        
            m=mm-NBx;
            n=nn-NBy;
            
            Epsr_PML(mm,nn)=epra*rect_2D(m,n,1,Tx,Ty,-f_pml*Tx/2,f_pml*Tx/2,  -Ty/2,Ty/2);
            Apsr_PML(mm,nn)=1/epra*rect_2D(m,n,1,Tx,Ty,-f_pml*Tx/2,f_pml*Tx/2,-Ty/2,Ty/2);
    end; % for nn
end; % for mm


Epsr_nonPML=Epsr_PML;
Apsr_nonPML=Apsr_PML;


for k=1:pml_N-1
    f_pml1=pml_fill_factor(k);
    f_pml2=pml_fill_factor(k+1);
    for mm=1:num_hx
        for nn=1:num_hy
            m=mm-NBx;
            n=nn-NBy;
            Epsr_PML(mm,nn)=Epsr_PML(mm,nn)+pml_epsr(k)*...
                (rect_2D(m,n,1,Tx,Ty,-f_pml2*Tx/2,f_pml2*Tx/2,-Ty/2,Ty/2) - rect_2D(m,n,1,Tx,Ty,-f_pml1*Tx/2,f_pml1*Tx/2,-Ty/2,Ty/2)  );
            Apsr_PML(mm,nn)=Apsr_PML(mm,nn)+1/pml_epsr(k)*...
                (rect_2D(m,n,1,Tx,Ty,-f_pml2*Tx/2,f_pml2*Tx/2,-Ty/2,Ty/2) - rect_2D(m,n,1,Tx,Ty,-f_pml1*Tx/2,f_pml1*Tx/2,-Ty/2,Ty/2)  );
        end;
    end;
end;

    f_pml1=pml_fill_factor(pml_N);
    f_pml2=1;
for mm=1:num_hx
        for nn=1:num_hy
            m=mm-NBx;
            n=nn-NBy;
            Epsr_PML(mm,nn)=Epsr_PML(mm,nn)+pml_epsr(pml_N)*...
                (rect_2D(m,n,1,Tx,Ty,-f_pml2*Tx/2,f_pml2*Tx/2,-Ty/2,Ty/2) - rect_2D(m,n,1,Tx,Ty,-f_pml1*Tx/2,f_pml1*Tx/2,-Ty/2,Ty/2)  );
            Apsr_PML(mm,nn)=Apsr_PML(mm,nn)+1/pml_epsr(pml_N)*...
                (rect_2D(m,n,1,Tx,Ty,-f_pml2*Tx/2,f_pml2*Tx/2,-Ty/2,Ty/2) - rect_2D(m,n,1,Tx,Ty,-f_pml1*Tx/2,f_pml1*Tx/2,-Ty/2,Ty/2)  );
        end;
 end;
 
 
 
%  
%  
% %%%  
% dx=Tx/num_hx;
% dy=Ty/num_hy;
% x=(-Tx/2+0.5*dx:dx:Tx/2-0.5*dx);
% y=(-Ty/2+0.5*dy:dy:Ty/2-0.5*dy);
% 
%
% Fx=2*pi/dx;
% Fy=2*pi/dy;
% dfx=Fx/num_hx;
% dfy=Fy/num_hy;
% fx=(-Fx/2+0.5*dfx:dfx:Fx/2-0.5*dfx);
% fy=(-Fy/2+0.5*dfy:dfy:Fy/2-0.5*dfy);
%     
% Gr_str=zeros(num_hx,num_hy);
% [xa ya]=meshgrid(x,y);
% for k=-2*nx:2*nx
%   for l=-2*ny:2*ny
%   Gr_str=Gr_str+Epsr_PML(k+NBx,l+NBy)*exp(j*(k*xa*2*pi/Tx+l*ya*2*pi/Ty));   
%   end;
% end;
% 
% figure(lnt);
% mesh(abs(Gr_str));
%  
%  
 
 
 
 
 
 