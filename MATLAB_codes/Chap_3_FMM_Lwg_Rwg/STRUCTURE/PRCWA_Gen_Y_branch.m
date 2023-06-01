% PRCWA_Gen_SPP_beaming.m
% In Grating_gen_Y_branch.m, the Foureir series coefficients of triangle
% grating structure with width of 0.5um and height of 1um are obtained.

% Layer structure ------------------------------------------
Nlay=30+2;

% material permittivity
epra=1;                       % air
eprm=(fun_Ag_nk(lambda))^2;   % gold
%eprm=2.5;

grating_thick=3*um;            % grating region

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



Grating_gen_Y_branch;        % Fourier series of grating structure
alpha_tm=1;
beta_tm=1;


%---------------------- e_h, a_h, g_h, b_h -----------------------
eps_xx=zeros(num_hx,num_hy,Nlay); eps_xy=zeros(num_hx,num_hy,Nlay); eps_xz=zeros(num_hx,num_hy,Nlay);   % permittivity
eps_yx=zeros(num_hx,num_hy,Nlay); eps_yy=zeros(num_hx,num_hy,Nlay); eps_yz=zeros(num_hx,num_hy,Nlay);
eps_zx=zeros(num_hx,num_hy,Nlay); eps_zy=zeros(num_hx,num_hy,Nlay); eps_zz=zeros(num_hx,num_hy,Nlay); 

aps_xx=zeros(num_hx,num_hy,Nlay); aps_xy=zeros(num_hx,num_hy,Nlay); aps_xz=zeros(num_hx,num_hy,Nlay);   % reciprocal permittivity 
aps_yx=zeros(num_hx,num_hy,Nlay); aps_yy=zeros(num_hx,num_hy,Nlay); aps_yz=zeros(num_hx,num_hy,Nlay);
aps_zx=zeros(num_hx,num_hy,Nlay); aps_zy=zeros(num_hx,num_hy,Nlay); aps_zz=zeros(num_hx,num_hy,Nlay); 

mu_xx=zeros(num_hx,num_hy,Nlay);  mu_xy=zeros(num_hx,num_hy,Nlay);  mu_xz=zeros(num_hx,num_hy,Nlay);    % permeability
mu_yx=zeros(num_hx,num_hy,Nlay);  mu_yy=zeros(num_hx,num_hy,Nlay);  mu_yz=zeros(num_hx,num_hy,Nlay);
mu_zx=zeros(num_hx,num_hy,Nlay);  mu_zy=zeros(num_hx,num_hy,Nlay);  mu_zz=zeros(num_hx,num_hy,Nlay); 

bu_xx=zeros(num_hx,num_hy,Nlay);  bu_xy=zeros(num_hx,num_hy,Nlay);  bu_xz=zeros(num_hx,num_hy,Nlay);    % reciprocal permeability
bu_yx=zeros(num_hx,num_hy,Nlay);  bu_yy=zeros(num_hx,num_hy,Nlay);  bu_yz=zeros(num_hx,num_hy,Nlay);
bu_zx=zeros(num_hx,num_hy,Nlay);  bu_zy=zeros(num_hx,num_hy,Nlay);  bu_zz=zeros(num_hx,num_hy,Nlay); 


for lnt=1:Nlay
    
       eps_xx(:,:,lnt)=eps_L(:,:,lnt);
       eps_yy(:,:,lnt)=eps_L(:,:,lnt);
       eps_zz(:,:,lnt)=eps_L(:,:,lnt);
       aps_xx(:,:,lnt)=aps_L(:,:,lnt);
       aps_yy(:,:,lnt)=aps_L(:,:,lnt);
       aps_zz(:,:,lnt)=aps_L(:,:,lnt);
       mu_xx(:,:,lnt)=  mu_L(:,:,lnt);
       mu_yy(:,:,lnt)=  mu_L(:,:,lnt);
       mu_zz(:,:,lnt)=  mu_L(:,:,lnt);
       bu_xx(:,:,lnt)=  bu_L(:,:,lnt);
       bu_yy(:,:,lnt)=  bu_L(:,:,lnt);
       bu_zz(:,:,lnt)=  bu_L(:,:,lnt);
               
end; % for lnt

%% permittivity profile testing        

% 2D
% xx=[-Tx/2:Tx*0.01:Tx/2]';
% Gr_str=zeros(length(xx),Nlay);    
% 
%   for lnt=1:Nlay
%      for k=-2*nx:2*nx
%          
%              Gr_str(:,lnt)=Gr_str(:,lnt)+eps_xx(k+NBx,1,lnt)*exp(j*(k*xx*2*pi/Tx)); 
%      end;
%   end;
%      figure(5);mesh(abs(Gr_str));

% 3D structure x-y crosssection
% xx=5*[-Tx/2:Tx*0.01:Tx/2];
% yy=5*[-Ty/2:Ty*0.01:Ty/2];
% 
%      Gr_str=zeros(length(xx),length(yy));    
%      [ya xa]=meshgrid(yy,xx);
%      
%      for k=-2*nx:2*nx
%          for l=-2*ny:2*ny
%      
%              Gr_str=Gr_str+eps_xx(k+NBx,l+NBy,29)*exp(j*(k*xa*2*pi/Tx+l*ya*2*pi/Ty)); 
%          
%          end;
%      end;
%      
%      figure(5);mesh(abs(Gr_str));
%%

