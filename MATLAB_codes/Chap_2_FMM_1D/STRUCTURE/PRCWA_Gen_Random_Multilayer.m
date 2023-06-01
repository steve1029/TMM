% PRCWA_Gen_SPP_beaming.m
% In Grating_gen_Y_branch.m, the Foureir series coefficients of triangle
% grating structure with width of 0.5um and height of 1um are obtained.

% Layer structure ------------------------------------------
Nlay=1+0;

% material permittivity
grating_thick=5*um;            % grating region

lay_thick=zeros(Nlay,1);      %   layer 
lay_thick(1)=0;                 %   left semi-infinite waveguide
for laynt=2:Nlay-1
  lay_thick(laynt)=grating_thick/(Nlay-2);  
end;
lay_thick(Nlay)=0;                 %   right semi-infinite waveguide

ac_thick=zeros(1,Nlay); 
for cnt1=1:Nlay
	 cnt2=1;
    while cnt2 <= cnt1
       ac_thick(cnt1)=ac_thick(cnt1)+lay_thick(cnt2);
       cnt2=cnt2+1;
    end;
end; % for cnt1


Grating_gen_Random_Multilayer;        % Permittivity & Permeability
alpha_tm=0;
beta_tm=0;


%---------------------- e_h, a_h, g_h, b_h -----------------------
eps_xx=zeros(1,Nlay); eps_xy=zeros(1,Nlay); eps_xz=zeros(1,Nlay);   % permittivity
eps_yx=zeros(1,Nlay); eps_yy=zeros(1,Nlay); eps_yz=zeros(1,Nlay);
eps_zx=zeros(1,Nlay); eps_zy=zeros(1,Nlay); eps_zz=zeros(1,Nlay); 

aps_xx=zeros(1,Nlay); aps_xy=zeros(1,Nlay); aps_xz=zeros(1,Nlay);   % reciprocal permittivity 
aps_yx=zeros(1,Nlay); aps_yy=zeros(1,Nlay); aps_yz=zeros(1,Nlay);
aps_zx=zeros(1,Nlay); aps_zy=zeros(1,Nlay); aps_zz=zeros(1,Nlay); 

mu_xx=zeros(1,Nlay);  mu_xy=zeros(1,Nlay);  mu_xz=zeros(1,Nlay);    % permeability
mu_yx=zeros(1,Nlay);  mu_yy=zeros(1,Nlay);  mu_yz=zeros(1,Nlay);
mu_zx=zeros(1,Nlay);  mu_zy=zeros(1,Nlay);  mu_zz=zeros(1,Nlay); 

bu_xx=zeros(1,Nlay);  bu_xy=zeros(1,Nlay);  bu_xz=zeros(1,Nlay);    % reciprocal permeability
bu_yx=zeros(1,Nlay);  bu_yy=zeros(1,Nlay);  bu_yz=zeros(1,Nlay);
bu_zx=zeros(1,Nlay);  bu_zy=zeros(1,Nlay);  bu_zz=zeros(1,Nlay); 


for lnt=1:Nlay    %diagonal anisotropic material multilayer
    
       eps_xx(1,lnt)=eps_L(1,lnt);
       eps_yy(1,lnt)=eps_L(1,lnt);
       eps_zz(1,lnt)=eps_L(1,lnt);
       aps_xx(1,lnt)=aps_L(1,lnt);
       aps_yy(1,lnt)=aps_L(1,lnt);
       aps_zz(1,lnt)=aps_L(1,lnt);
       mu_xx(1,lnt)=  mu_L(1,lnt);
       mu_yy(1,lnt)=  mu_L(1,lnt);
       mu_zz(1,lnt)=  mu_L(1,lnt);
       bu_xx(1,lnt)=  bu_L(1,lnt);
       bu_yy(1,lnt)=  bu_L(1,lnt);
       bu_zz(1,lnt)=  bu_L(1,lnt);
               
end; % for lnt


