% SPP beaming structure
 
  Wx=100*nm;
  
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
   
   
   for lnt=1:Nlay
        

    if lnt == 1
        %----------------------------------------------- SLIT REGION
                                         
        eps=epra*rect_2D_mesh(m,n,1,Tx,Ty,-Wx/2,Wx/2,-Ty/2,Ty/2)+eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-Wx/2,Wx/2,-Ty/2,Ty/2));
        aps=1/epra*rect_2D_mesh(m,n,1,Tx,Ty,-Wx/2,Wx/2,-Ty/2,Ty/2)+1/eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-Wx/2,Wx/2,-Ty/2,Ty/2));
        mu=mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
        bu=1/mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);  
        
        eps_L(:,:,lnt)=eps;
        aps_L(:,:,lnt)=aps;
        mu_L(:,:,lnt)=mu;
        bu_L(:,:,lnt)=bu;
    
    end; % if lnt
        
     if lnt == 2
        %---------------------------------------------- Dielectric grating
        Wxg=Wx+100*nano;
        KK=8;
        facR=0.5;
        facL=0.5;
        eprgR=1.72^2;
        eprgL=eprgR;
        TgR=380*nm;
        TgL=TgR;
%         TgR=305*nm;
%         TgL=505*nm;
        
        
      % permittivity   
        eps=epra*rect_2D_mesh(m,n,1,Tx,Ty,-Wxg/2,Wxg/2,-Ty/2,Ty/2);
        aps=1/epra*rect_2D_mesh(m,n,1,Tx,Ty,-Wxg/2,Wxg/2,-Ty/2,Ty/2);
                 
       for kk=0:KK
        eps= eps+( eprgR*rect_2D_mesh(m,n,1,Tx,Ty,Wxg/2+kk*TgR,    Wxg/2+(kk+facR)*TgR,-Ty/2,Ty/2) + eprgL*rect_2D_mesh(m,n,1,Tx,Ty,-Wxg/2-(kk+facL)*TgL,-Wxg/2-kk*TgL,-Ty/2,Ty/2)  ) ...
            +    ( epra *rect_2D_mesh(m,n,1,Tx,Ty,Wxg/2+(kk+1)*TgR,Wxg/2+(kk+facR)*TgR,-Ty/2,Ty/2) + epra *rect_2D_mesh(m,n,1,Tx,Ty,-Wxg/2-(kk+1)*TgL,   -Wxg/2-(kk+facL)*TgL,-Ty/2,Ty/2)  ) ;
              
        aps=  aps+( 1/eprgR*rect_2D_mesh(m,n,1,Tx,Ty,Wxg/2+kk*TgR,    Wxg/2+(kk+facR)*TgR,-Ty/2,Ty/2) + 1/eprgL*rect_2D_mesh(m,n,1,Tx,Ty,-Wxg/2-(kk+facL)*TgL,-Wxg/2-kk*TgL,-Ty/2,Ty/2)  ) ...
            +     ( 1/epra *rect_2D_mesh(m,n,1,Tx,Ty,Wxg/2+(kk+1)*TgR,Wxg/2+(kk+facR)*TgR,-Ty/2,Ty/2) + 1/epra *rect_2D_mesh(m,n,1,Tx,Ty,-Wxg/2-(kk+1)*TgL,   -Wxg/2-(kk+facL)*TgL,-Ty/2,Ty/2)  ) ;
        
        end;
        
         eps= eps+  epra*( rect_2D_mesh(m,n,1,Tx,Ty,Wxg/2+(KK+1)*TgR,Tx/2-pml_width,-Ty/2,Ty/2) + rect_2D_mesh(m,n,1,Tx,Ty,-Wxg/2-(KK+1)*TgL,-Tx/2+pml_width,-Ty/2,Ty/2)  ) ;
         aps= aps+1/epra*( rect_2D_mesh(m,n,1,Tx,Ty,Wxg/2+(KK+1)*TgR,Tx/2-pml_width,-Ty/2,Ty/2) + rect_2D_mesh(m,n,1,Tx,Ty,-Wxg/2-(KK+1)*TgL,-Tx/2+pml_width,-Ty/2,Ty/2)  ) ;

        mu=mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
        bu=1/mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);  
    
        eps=eps+ Epsr_PML-Epsr_nonPML;
        aps=aps+ Apsr_PML-Apsr_nonPML;
       
        eps_L(:,:,lnt)=eps;
        aps_L(:,:,lnt)=aps;
        mu_L(:,:,lnt)=mu;
        bu_L(:,:,lnt)=bu;
        
    end; % if lnt
   
    if lnt == 3
       %--------------------------------------------------- freespace+PML
        eps=Epsr_PML;
        aps=Apsr_PML;
        mu=mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
        bu=1/mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);  
        
        eps_L(:,:,lnt)=eps;
        aps_L(:,:,lnt)=aps;
        mu_L(:,:,lnt)=mu;
        bu_L(:,:,lnt)=bu;
        
              
   end; % if lnt
    
end; % for lnt
   