% diagonal anisotropic material

  str_tensor=struct('eps_xx',eps_xx(:,:,laynt), 'eps_yy',eps_yy(:,:,laynt),'eps_zz',eps_zz(:,:,laynt), ...
                    'aps_xx',aps_xx(:,:,laynt),'aps_yy',aps_yy(:,:,laynt), 'aps_zz',aps_zz(:,:,laynt), ...
                    'mu_xx',  mu_xx(:,:,laynt),'mu_yy',  mu_yy(:,:,laynt),'mu_zz',  mu_zz(:,:,laynt), ...
                    'bu_xx',  bu_xx(:,:,laynt),'bu_yy',  bu_yy(:,:,laynt),'bu_zz',  bu_zz(:,:,laynt));
  
  [T_tempa, R_tempa, T_tempb, R_tempb, Laycofa, Laycofb,pfEx,pfEy,pfEz,pfHx,pfHy,pfHz,pevalue,mfEx,mfEy,mfEz,mfHx,mfHy,mfHz,mevalue]=FMM_single_block(lay_thick(laynt),str_tensor,alpha_tm,beta_tm);

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
       