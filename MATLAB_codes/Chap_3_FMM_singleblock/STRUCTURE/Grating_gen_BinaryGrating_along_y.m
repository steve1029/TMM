% Binary grating

   Wy=0.5*um;
   Height=1*um;

  [m,n]=meshgrid([1:1:num_hx]-NBx,[1:1:num_hy]-NBy);
   m=m';
   n=n';
   
   
   Wy_L=Wy;
   
   eps_L=zeros(num_hx,num_hy,Nlay);
   aps_L=zeros(num_hx,num_hy,Nlay);
   mu_L =zeros(num_hx,num_hy,Nlay);
   bu_L =zeros(num_hx,num_hy,Nlay);
   
   for lnt=1:Nlay
   
   wy=Wy_L(lnt);
       
   eps=eprb  *rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-wy/2,wy/2)   + eprm*( rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-wy/2,wy/2));
   aps=1/eprb*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-wy/2,wy/2) + 1/eprm*( rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-wy/2,wy/2));  
   mu=mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
   bu=1/mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);  

   eps_L(:,:,lnt)=eps;
   aps_L(:,:,lnt)=aps;
   mu_L(:,:,lnt)=mu;
   bu_L(:,:,lnt)=bu;
   
   end; % 
   
   
   