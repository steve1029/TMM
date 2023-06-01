% Grating_gen_PCWG

xc=[-7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7]*a;

r=0.2*a; % = 6*dx
for m=1:15
    for k=1:length(sx)
     for l=1:length(sz)
        if r^2-(sx(k)-xc(m))^2-(sz(l))^2 > 0  &  xc(m) ~=0
        Str(k,l)=1;
        end;
    end;
end;
end;

figure(34);
imagesc(sz/a,sx/a,abs(Str)); 

[m,n]=meshgrid([1:1:num_hx]-NBx,[1:1:num_hy]-NBy);
 m=m';
 n=n';

Aeps_L=zeros(num_hx,num_hy,Nlay);
Aaps_L=zeros(num_hx,num_hy,Nlay);
Amu_L =zeros(num_hx,num_hy,Nlay );
Abu_L =zeros(num_hx,num_hy,Nlay );



Gr_str=zeros(length(sx),length(sz)); % grating structure

[m,n]=meshgrid([1:1:num_hx]-NBx,[1:1:num_hy]-NBy);
 m=m';
 n=n';
        
for lnt=1:Nlay
                
                eps=zeros(num_hx,num_hy);
                aps=zeros(num_hx,num_hy);
                mu =zeros(num_hx,num_hy);
                bu =zeros(num_hx,num_hy);

                for k=1:length(sx)  % for x-axis 
                                    
                    if Str(k,lnt) == 1
                    eps=eps+eprm*rect_2D_mesh(m,n,1,Tx,Ty,sx(k)-dx/2,sx(k)+dx/2,-Ty/2,Ty/2);
                    aps=aps+1/eprm*rect_2D_mesh(m,n,1,Tx,Ty,sx(k)-dx/2,sx(k)+dx/2,-Ty/2,Ty/2);  
                    else
                    eps=eps+epra*rect_2D_mesh(m,n,1,Tx,Ty,sx(k)-dx/2,sx(k)+dx/2,-Ty/2,Ty/2);
                    aps=aps+1/epra*rect_2D_mesh(m,n,1,Tx,Ty,sx(k)-dx/2,sx(k)+dx/2,-Ty/2,Ty/2);
                    end; % if Str
                    
                end;
                
                mu=mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
                bu=1/mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2); 


     
        Aeps_L(:,:,lnt)=eps;
        Aaps_L(:,:,lnt)=aps;
        Amu_L(:,:,lnt)=mu;
        Abu_L(:,:,lnt)=bu;

        
     for k=-2*nx:2*nx
         for l=-2*ny:2*ny
     
             Gr_str(:,lnt)=Gr_str(:,lnt)+eps(k+NBx,l+NBy)*exp(j*(k*sx'*2*pi/Tx)); 
         
         end;
     end;
        
end; % for lnt


figure(31);mesh(abs(Gr_str)); axis equal
        

