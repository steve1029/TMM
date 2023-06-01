% SPP beaming structure

Wx_A=400*nm;
Wx_B=100*nm;

[m,n]=meshgrid((1:1:num_hx)-NBx,(1:1:num_hy)-NBy);
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


    end; % if lnt

    if lnt == 2
        %----------------------------------------------- SLIT A REGION

        eps=epra*rect_2D_mesh(m,n,1,Tx,Ty,-Wx_A/2,Wx_A/2,-Ty/2,Ty/2)+eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-Wx_A/2,Wx_A/2,-Ty/2,Ty/2));
        aps=1/epra*rect_2D_mesh(m,n,1,Tx,Ty,-Wx_A/2,Wx_A/2,-Ty/2,Ty/2)+1/eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-Wx_A/2,Wx_A/2,-Ty/2,Ty/2));
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

    end; % if lnt

    if lnt == 3
        %----------------------------------------------- SLIT B REGION

        eps=epra*rect_2D_mesh(m,n,1,Tx,Ty,-Wx_B/2,Wx_B/2,-Ty/2,Ty/2)+eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-Wx_B/2,Wx_B/2,-Ty/2,Ty/2));
        aps=1/epra*rect_2D_mesh(m,n,1,Tx,Ty,-Wx_B/2,Wx_B/2,-Ty/2,Ty/2)+1/eprm*(rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2)-rect_2D_mesh(m,n,1,Tx,Ty,-Wx_B/2,Wx_B/2,-Ty/2,Ty/2));
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

    end; % if lnt


end; % for lnt
