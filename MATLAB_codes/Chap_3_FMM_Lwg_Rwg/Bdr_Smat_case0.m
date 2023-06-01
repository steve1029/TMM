% N-block interconnection
% free space - grating - free space       
 
%% part 5 Field visualization & data save

K=zeros(2*L,2*L); % KII matrix
for k=1:NBx
   for l=1:NBy
od_ind1=(k-1)*NBy+l;      
K11(od_ind1,od_ind1)=(kx_vc(k)*ky_vc(l))/kz_vc(k,l);
K12(od_ind1,od_ind1)=(kz_vc(k,l)^2+kx_vc(k)^2)/kz_vc(k,l);
K21(od_ind1,od_ind1)=-(ky_vc(l)^2+kz_vc(k,l)^2)/kz_vc(k,l);
K22(od_ind1,od_ind1)=-( ky_vc(l)*kx_vc(k) )/kz_vc(k,l);

   end;
end;

K(1:L, 1:L)=K11;
K(1:L, L+1:2*L)=K12;
K(L+1:2*L, 1:L)=K21;
K(L+1:2*L,L+1:2*L )=K22;

Wh=Iden;
Vh=K/(w0*mu0);

switch direct_

    case 1   % left-to-right
        
        Field_visualization_3D_xz_case0_Lfree_Rfree_leftright;

    case 2   % right-to-left
        
        Field_visualization_3D_xz_case0_Lfree_Rfree_rightleft;

end;

