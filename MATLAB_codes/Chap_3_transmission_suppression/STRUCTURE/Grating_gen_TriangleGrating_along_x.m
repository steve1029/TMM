% triangle grating
% Nlay

   Wx=0.5*um;
   Height=1*um;

 [m,n]=meshgrid([1:1:num_hx]-NBx,[1:1:num_hy]-NBy);
   m=m';
   n=n';
   
   
   Wx_L=[Wx:-Wx/(Nlay-1):0];
   
   eps_L=zeros(num_hx,num_hy,Nlay);
   aps_L=zeros(num_hx,num_hy,Nlay);
   mu_L =zeros(num_hx,num_hy,Nlay);
   bu_L =zeros(num_hx,num_hy,Nlay);
   
   for lnt=1:Nlay
   
   wx=Wx_L(lnt);
       
   eps=eprb  *rect_2D_mesh(m,n,1,Tx,Ty,-wx/2,wx/2,-Ty/2,Ty/2)   + eprm*( rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-wx/2,wx/2,-Ty/2,Ty/2));
   aps=1/eprb*rect_2D_mesh(m,n,1,Tx,Ty,-wx/2,wx/2,-Ty/2,Ty/2) + 1/eprm*( rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-wx/2,wx/2,-Ty/2,Ty/2));  
   mu=mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
   bu=1/mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);  

   eps_L(:,:,lnt)=eps;
   aps_L(:,:,lnt)=aps;
   mu_L(:,:,lnt)=mu;
   bu_L(:,:,lnt)=bu;
   
   end; % 
   
   
   