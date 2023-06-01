% SPP beaming structure
 
  Wx=Tx/6;  % waveguide width
  Sx=Tx/4;  % gap between waveguides
  
  [m,n]=meshgrid([1:1:num_hx]-NBx,[1:1:num_hy]-NBy);
        m=m';
        n=n';
         
   eps_L=zeros(num_hx,num_hy,Nlay);
   aps_L=zeros(num_hx,num_hy,Nlay);
   mu_L =zeros(num_hx,num_hy,Nlay);
   bu_L =zeros(num_hx,num_hy,Nlay);
   
   eps=zeros(num_hx,num_hy,Nlay);
   aps=zeros(num_hx,num_hy,Nlay);
   mu =zeros(num_hx,num_hy,Nlay);
   bu =zeros(num_hx,num_hy,Nlay);
   
   % left semi-infinite waveguide
   lnt = 1;
                                      
        eps=epra*rect_2D_mesh(m,n,1,Tx,Ty,-Wx/2,Wx/2,-Ty/2,Ty/2)+eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-Wx/2,Wx/2,-Ty/2,Ty/2));
        aps=1/epra*rect_2D_mesh(m,n,1,Tx,Ty,-Wx/2,Wx/2,-Ty/2,Ty/2)+1/eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-Wx/2,Wx/2,-Ty/2,Ty/2));
        mu=mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
        bu=1/mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);  
        
        eps_L(:,:,lnt)=eps;
        aps_L(:,:,lnt)=aps;
        mu_L(:,:,lnt)=mu;
        bu_L(:,:,lnt)=bu;
        
   % right semi-infinite waveguide 
  lnt=Nlay;  
    
        eps=eprm*rect_2D_mesh(m,n,1,Tx,Ty,-Sx/2,Sx/2,-Ty/2,Ty/2)...
            +epra*(rect_2D_mesh(m,n,1,Tx,Ty,Sx/2,Sx/2+Wx,-Ty/2,Ty/2) + rect_2D_mesh(m,n,1,Tx,Ty,-Sx/2,-Sx/2-Wx,-Ty/2,Ty/2))... 
            +eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-Sx/2-Wx,Sx/2+Wx,-Ty/2,Ty/2));
      
        aps=1/eprm*rect_2D_mesh(m,n,1,Tx,Ty,-Sx/2,Sx/2,-Ty/2,Ty/2)...
            +1/epra*(rect_2D_mesh(m,n,1,Tx,Ty,Sx/2,Sx/2+Wx,-Ty/2,Ty/2) + rect_2D_mesh(m,n,1,Tx,Ty,-Sx/2,-Sx/2-Wx,-Ty/2,Ty/2))... 
            +1/eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-Sx/2-Wx,Sx/2+Wx,-Ty/2,Ty/2));
        
        mu=mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
        bu=1/mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);  
        
        eps_L(:,:,lnt)=eps;
        aps_L(:,:,lnt)=aps;
        mu_L(:,:,lnt)=mu;
        bu_L(:,:,lnt)=bu;
   

   for lnt=2:Nlay-1
        
       if ac_thick(lnt) < grating_thick*Wx/(Sx+Wx)
        
        gwth=(Sx+Wx)/grating_thick*ac_thick(lnt)+Wx;
        
        eps=epra*rect_2D_mesh(m,n,1,Tx,Ty,-gwth/2,gwth/2,-Ty/2,Ty/2)+eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-gwth/2,gwth/2,-Ty/2,Ty/2));
        aps=1/epra*rect_2D_mesh(m,n,1,Tx,Ty,-gwth/2,gwth/2,-Ty/2,Ty/2)+1/eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-gwth/2,gwth/2,-Ty/2,Ty/2));
        mu=mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
        bu=1/mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);  
        
        eps_L(:,:,lnt)=eps;
        aps_L(:,:,lnt)=aps;
        mu_L(:,:,lnt)=mu;
        bu_L(:,:,lnt)=bu;
         
       else
           
        gsth=(Sx+Wx)/grating_thick*ac_thick(lnt)-Wx;
           
        eps=eprm*rect_2D_mesh(m,n,1,Tx,Ty,-gsth/2,gsth/2,-Ty/2,Ty/2)...
            +epra*(rect_2D_mesh(m,n,1,Tx,Ty,gsth/2,gsth/2+Wx,-Ty/2,Ty/2) + rect_2D_mesh(m,n,1,Tx,Ty,-gsth/2,-gsth/2-Wx,-Ty/2,Ty/2))... 
            +eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-gsth/2-Wx,gsth/2+Wx,-Ty/2,Ty/2));
      
        aps=1/eprm*rect_2D_mesh(m,n,1,Tx,Ty,-gsth/2,gsth/2,-Ty/2,Ty/2)...
            +1/epra*(rect_2D_mesh(m,n,1,Tx,Ty,gsth/2,gsth/2+Wx,-Ty/2,Ty/2) + rect_2D_mesh(m,n,1,Tx,Ty,-gsth/2,-gsth/2-Wx,-Ty/2,Ty/2))... 
            +1/eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-gsth/2-Wx,gsth/2+Wx,-Ty/2,Ty/2));
        
        mu=mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
        bu=1/mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);  
        
        eps_L(:,:,lnt)=eps;
        aps_L(:,:,lnt)=aps;
        mu_L(:,:,lnt)=mu;
        bu_L(:,:,lnt)=bu;
           
       end;
    

end; % for lnt
   