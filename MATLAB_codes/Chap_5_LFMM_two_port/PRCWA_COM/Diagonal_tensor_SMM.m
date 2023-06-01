% diagonal anisotropic material

eps_xx=Aeps_xx; eps_xy=Aeps_xy; eps_xz=Aeps_xz;
eps_yx=Aeps_yx; eps_yy=Aeps_yy; eps_yz=Aeps_yz;
eps_zx=Aeps_zx; eps_zy=Aeps_zy; eps_zz=Aeps_zz;

aps_xx=Aaps_xx; aps_xy=Aaps_xy; aps_xz=Aaps_xz;
aps_yx=Aaps_yx; aps_yy=Aaps_yy; aps_yz=Aaps_yz;
aps_zx=Aaps_zx; aps_zy=Aaps_zy; aps_zz=Aaps_zz;

mu_xx=Amu_xx; mu_xy=Amu_xy; mu_xz=Amu_xz;
mu_yx=Amu_yx; mu_yy=Amu_yy; mu_yz=Amu_yz;
mu_zx=Amu_zx; mu_zy=Amu_zy; mu_zz=Amu_zz;

bu_xx=Abu_xx; bu_xy=Abu_xy; bu_xz=Abu_xz;
bu_yx=Abu_yx; bu_yy=Abu_yy; bu_yz=Abu_yz;
bu_zx=Abu_zx; bu_zy=Abu_zy; bu_zz=Abu_zz;


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
 
 