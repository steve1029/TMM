% diagonal anisotropic material

for laynt=1:Nlay

         % diagonal anisotropic material
        str_tensor=struct('eps_xx',eps_xx(:,:,laynt), 'eps_yy',eps_yy(:,:,laynt),'eps_zz',eps_zz(:,:,laynt), ...
                    'aps_xx',aps_xx(:,:,laynt),'aps_yy',aps_yy(:,:,laynt), 'aps_zz',aps_zz(:,:,laynt), ...
                    'mu_xx',  mu_xx(:,:,laynt),'mu_yy',  mu_yy(:,:,laynt),'mu_zz',  mu_zz(:,:,laynt), ...
                    'bu_xx',  bu_xx(:,:,laynt),'bu_yy',  bu_yy(:,:,laynt),'bu_zz',  bu_zz(:,:,laynt));
  
        [T_tempa, R_tempa, T_tempb, R_tempb, Laycofa, Laycofb,pfEx,pfEy,pfEz,pfHx,pfHy,pfHz,pevalue,mfEx,mfEy,mfEz,mfHx,mfHy,mfHz,mevalue]=FMM_single_block(lay_thick(laynt),str_tensor,alpha_tm,beta_tm);

    
    Ta(:,:,laynt)=T_tempa;
    Ra(:,:,laynt)=R_tempa;
    Tb(:,:,laynt)=T_tempb;
    Rb(:,:,laynt)=R_tempb;
    Ca(:,:,laynt)=Laycofa;
    Cb(:,:,laynt)=Laycofb;
    
    Pf_Ex(:,:,laynt)=pfEx;
    Pf_Ey(:,:,laynt)=pfEy;
    Pf_Ez(:,:,laynt)=pfEz;
 
    Pf_Hx(:,:,laynt)=pfHx;
    Pf_Hy(:,:,laynt)=pfHy;
    Pf_Hz(:,:,laynt)=pfHz;

    Mf_Ex(:,:,laynt)=mfEx;
    Mf_Ey(:,:,laynt)=mfEy;
    Mf_Ez(:,:,laynt)=mfEz;
 
    Mf_Hx(:,:,laynt)=mfHx;
    Mf_Hy(:,:,laynt)=mfHy;
    Mf_Hz(:,:,laynt)=mfHz;
    
    Peigvalue(1,:,laynt)=pevalue;
    Meigvalue(1,:,laynt)=mevalue;
       
 end;  % for laynt
 
 
% Boundary S-matrix computation of Input Output

% Bdr_SMat_infr_outfr : free-space boundary S-matrix 
% Bdr_SMat_wg_tensor  : waveguide boundary S-matrix


PRCWA_Gen_inout_K;                          %   input-output region 
[Lfree_Tf,Lfree_Rb,Lfree_Tb,Lfree_Rf,Rfree_Tf,Rfree_Rb,Rfree_Tb,Rfree_Rf]=Bdr_SMat_infr_outfr();

laynt=1;  %  left side semi-infinite waveguide 
% diagonal anisotropic material
 str_tensor=struct('eps_xx',eps_xx(:,:,laynt), 'eps_yy',eps_yy(:,:,laynt),'eps_zz',eps_zz(:,:,laynt), ...
                    'aps_xx',aps_xx(:,:,laynt),'aps_yy',aps_yy(:,:,laynt), 'aps_zz',aps_zz(:,:,laynt), ...
                    'mu_xx',  mu_xx(:,:,laynt),'mu_yy',  mu_yy(:,:,laynt),'mu_zz',  mu_zz(:,:,laynt), ...
                    'bu_xx',  bu_xx(:,:,laynt),'bu_yy', bu_yy(:,:,laynt),'bu_zz',  bu_zz(:,:,laynt));

[Lwg_Tf1,Lwg_Rb1,Lwg_Tb1,Lwg_Rf1,Rwg_Tf1,Rwg_Rb1,Rwg_Tb1,Rwg_Rf1]=Bdr_SMat_wg(str_tensor,alpha_tm,beta_tm);


laynt=Nlay;   % right side semi-infinite waveguide
% diagonal anisotropic material
 str_tensor=struct('eps_xx',eps_xx(:,:,laynt), 'eps_yy',eps_yy(:,:,laynt),'eps_zz',eps_zz(:,:,laynt), ...
                    'aps_xx',aps_xx(:,:,laynt),'aps_yy',aps_yy(:,:,laynt), 'aps_zz',aps_zz(:,:,laynt), ...
                    'mu_xx',  mu_xx(:,:,laynt),'mu_yy',  mu_yy(:,:,laynt),'mu_zz',  mu_zz(:,:,laynt), ...
                    'bu_xx',  bu_xx(:,:,laynt),'bu_yy',  bu_yy(:,:,laynt),'bu_zz',  bu_zz(:,:,laynt));

[Lwg_Tf2,Lwg_Rb2,Lwg_Tb2,Lwg_Rf2,Rwg_Tf2,Rwg_Rb2,Rwg_Tb2,Rwg_Rf2]=Bdr_SMat_wg(str_tensor,alpha_tm,beta_tm);


