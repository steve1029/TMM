% full off-diagonal anisotropic material

str_tensor=struct('eps_xx',eps_xx(:,:,laynt),'eps_xy',eps_xy(:,:,laynt),'eps_xz',eps_xz(:,:,laynt), ...
                         'eps_yx',eps_yx(:,:,laynt),'eps_yy',eps_yy(:,:,laynt),'eps_yz',eps_yz(:,:,laynt), ...
                         'eps_zx',eps_zx(:,:,laynt),'eps_zy',eps_zy(:,:,laynt),'eps_zz',eps_zz(:,:,laynt), ...
                         'aps_xx',aps_xx(:,:,laynt),'aps_xy',aps_xy(:,:,laynt),'aps_xz',aps_xz(:,:,laynt), ...
                         'aps_yx',aps_yx(:,:,laynt),'aps_yy',aps_yy(:,:,laynt),'aps_yz',aps_yz(:,:,laynt), ...
                         'aps_zx',aps_zx(:,:,laynt),'aps_zy',aps_zy(:,:,laynt),'aps_zz',aps_zz(:,:,laynt), ...
                         'mu_xx',  mu_xx(:,:,laynt),'mu_xy',  mu_xy(:,:,laynt),'mu_xz',  mu_xz(:,:,laynt), ...
                         'mu_yx',  mu_yx(:,:,laynt),'mu_yy',  mu_yy(:,:,laynt),'mu_yz',  mu_yz(:,:,laynt), ...
                         'mu_zx',  mu_zx(:,:,laynt),'mu_zy',  mu_zy(:,:,laynt),'mu_zz',  mu_zz(:,:,laynt), ...
                         'bu_xx',  bu_xx(:,:,laynt),'bu_xy',  bu_xy(:,:,laynt),'bu_xz',  bu_xz(:,:,laynt), ...
                         'bu_yx',  bu_yx(:,:,laynt),'bu_yy',  bu_yy(:,:,laynt),'bu_yz',  bu_yz(:,:,laynt), ...
                         'bu_zx',  bu_zx(:,:,laynt),'bu_zy',  bu_zy(:,:,laynt),'bu_zz',  bu_zz(:,:,laynt) );
                     
   [T_tempa, R_tempa, T_tempb, R_tempb, Laycofa, Laycofb,pfEx,pfEy,pfEz,pfHx,pfHy,pfHz,pevalue,mfEx,mfEy,mfEz,mfHx,mfHy,mfHz,mevalue]=FMM_single_block_tensor(lay_thick(laynt),str_tensor,alpha_tm,beta_tm);
   

    TTa=T_tempa;
    RRa=R_tempa;
    TTb=T_tempb;
    RRb=R_tempb;
    Ca=Laycofa;
    Cb=Laycofb;
    
    Pf_Ex=pfEx;
    Pf_Ey=pfEy;
    Pf_Ez=pfEz;
 
    Pf_Hx=pfHx;
    Pf_Hy=pfHy;
    Pf_Hz=pfHz;

    Mf_Ex=mfEx;
    Mf_Ey=mfEy;
    Mf_Ez=mfEz;
 
    Mf_Hx=mfHx;
    Mf_Hy=mfHy;
    Mf_Hz=mfHz;
    
    Peigvalue=pevalue;
    Meigvalue=mevalue;