% SPP beaming structure

Wx=100*nm;

[m,n]=meshgrid([1:1:num_hx]-NBx,[1:1:num_hy]-NBy);
m=m';
n=n';


% % no or old pml
% eps_L=zeros(num_hx,num_hy,Nlay);
% aps_L=zeros(num_hx,num_hy,Nlay);
% mu_L=zeros(num_hx,num_hy,Nlay);
% bu_L=zeros(num_hx,num_hy,Nlay);

% new pml
eps_L_xx=zeros(num_hx,num_hy,Nlay);
eps_L_yy=zeros(num_hx,num_hy,Nlay);
eps_L_zz=zeros(num_hx,num_hy,Nlay);
aps_L_xx=zeros(num_hx,num_hy,Nlay);
aps_L_yy=zeros(num_hx,num_hy,Nlay);
aps_L_zz=zeros(num_hx,num_hy,Nlay);
mu_L_xx =zeros(num_hx,num_hy,Nlay);
mu_L_yy =zeros(num_hx,num_hy,Nlay);
mu_L_zz =zeros(num_hx,num_hy,Nlay);
bu_L_xx =zeros(num_hx,num_hy,Nlay);
bu_L_yy =zeros(num_hx,num_hy,Nlay);
bu_L_zz =zeros(num_hx,num_hy,Nlay);

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

%         % no or old pml
%         eps_L(:,:,Int)=eps;
%         aps_L(:,:,Int)=aps;
%         mu_L(:,:,Int)=mu;
%         bu_L(:,:,Int)=bu;
        
        % new pml
        eps_L_xx(:,:,lnt)=eps;
        eps_L_yy(:,:,lnt)=eps;
        eps_L_zz(:,:,lnt)=eps;
        aps_L_xx(:,:,lnt)=aps;
        aps_L_yy(:,:,lnt)=aps;
        aps_L_zz(:,:,lnt)=aps;
        mu_L_xx(:,:,lnt)=mu;
        mu_L_yy(:,:,lnt)=mu;
        mu_L_zz(:,:,lnt)=mu;
        bu_L_xx(:,:,lnt)=bu;
        bu_L_yy(:,:,lnt)=bu;
        bu_L_zz(:,:,lnt)=bu;
figure(lnt);imagesc(abs(bu_L_zz(:,:,lnt)));title([num2str(lnt) ' / ' num2str(Nlay)]);
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

%         % no pml
%         eps_L(:,:,Int)=eps;
%         aps_L(:,:,Int)=aps;
%         mu_L(:,:,Int)=mu;
%         bu_L(:,:,Int)=bu;

        
%         % old pml
%         eps=eps+ Epsr_PML-Epsr_nonPML;
%         aps=aps+ Apsr_PML-Apsr_nonPML;
%         eps_L(:,:,Int)=eps;
%         aps_L(:,:,Int)=aps;
%         mu_L(:,:,Int)=mu;
%         bu_L(:,:,Int)=bu;

        % new pml
        eps_L_xx(:,:,lnt)=eps + Epsr_PML_xx;
        eps_L_yy(:,:,lnt)=eps + Epsr_PML_yy;
        eps_L_zz(:,:,lnt)=eps + Epsr_PML_zz;
        aps_L_xx(:,:,lnt)=aps + Apsr_PML_xx;
        aps_L_yy(:,:,lnt)=aps + Apsr_PML_yy;
        aps_L_zz(:,:,lnt)=aps + Apsr_PML_zz;
        mu_L_xx (:,:,lnt)=mu  + Mpsr_PML_xx;
        mu_L_yy (:,:,lnt)=mu  + Mpsr_PML_yy;
        mu_L_zz (:,:,lnt)=mu  + Mpsr_PML_zz;
        bu_L_xx (:,:,lnt)=bu  + Bpsr_PML_xx;
        bu_L_yy (:,:,lnt)=bu  + Bpsr_PML_yy;
        bu_L_zz (:,:,lnt)=bu  + Bpsr_PML_zz;
figure(lnt);imagesc(abs(bu_L_zz(:,:,lnt)));title([num2str(lnt) ' / ' num2str(Nlay)]);
    end; % if lnt

    if lnt == 3
        %--------------------------------------------------- freespace+PML
        %         eps=eps0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
        %         eps=eps+Epsr_PML;
        %         aps=Apsr_PML;
        %         mu=mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
        %         bu=1/mur0*rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);
        %
        % %         eps_L(:,:,lnt)=eps;
        % %         aps_L(:,:,lnt)=aps;

        %         % no PML
        %         eps_L(:,:,lnt)=mu;
        %         aps_L(:,:,lnt)=bu;
        %         mu_L(:,:,lnt)=mu;
        %         bu_L(:,:,lnt)=bu;

        % new PML
        eps_L_xx(:,:,lnt)=  epra*rect_region0 + Epsr_PML_xx;
        eps_L_yy(:,:,lnt)=  epra*rect_region0 + Epsr_PML_yy;
        eps_L_zz(:,:,lnt)=  epra*rect_region0 + Epsr_PML_zz;
        aps_L_xx(:,:,lnt)=1/epra*rect_region0 + Apsr_PML_xx;
        aps_L_yy(:,:,lnt)=1/epra*rect_region0 + Apsr_PML_yy;
        aps_L_zz(:,:,lnt)=1/epra*rect_region0 + Apsr_PML_zz;
        mu_L_xx (:,:,lnt)=  mur0*rect_region0 + Mpsr_PML_xx;
        mu_L_yy (:,:,lnt)=  mur0*rect_region0 + Mpsr_PML_yy;
        mu_L_zz (:,:,lnt)=  mur0*rect_region0 + Mpsr_PML_zz;
        bu_L_xx (:,:,lnt)=1/mur0*rect_region0 + Bpsr_PML_xx;
        bu_L_yy (:,:,lnt)=1/mur0*rect_region0 + Bpsr_PML_yy;
        bu_L_zz (:,:,lnt)=1/mur0*rect_region0 + Bpsr_PML_zz;

figure(lnt);imagesc(abs(bu_L_zz(:,:,lnt)));title([num2str(lnt) ' / ' num2str(Nlay)]);
    end; % if lnt

end; % for lnt
