% single binary grating structure analysis

% Layer structure
Nlay=1;    %  number of layers

% material permittivity
epra=1+0.00000001*i;                       % air
eprb=(1.4+0.1*i)^2;
eprm=(fun_Ag_nk(lambda))^2;   % gold
%eprm=(1.7+0.000000001*i)^2;


%Open Grating_gen_BinaryGrating_along_x.m and
%Grating_gen_BinaryGrating_along_y.m 
% % Let us learn the use of rect_2D_mesh.m
%rect_2D_mesh(m,n,1,Tx,Ty,-wx/2,wx/2,-Ty/2,Ty/2) produce the Fourier series
%of the x-directionally varying periodic binary function with the period of
%Tx and the width of wx. This function actually calculates the 1D Fourier 
%series of the 1D periodic rectangle fuction.

%rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-wy/2,wy/2) produce the Fourier series
%of the y-directionally varying periodic binary function with the period of
%Ty and the width of wy. This function actually calculates the 1D Fourier 
%series of the 1D periodic rectangle fuction.

%rect_2D_mesh(m,n,1,Tx,Ty,-wx/2,wx/2,-wy/2,wy/2) produce the Fourier series
%of the x-directionally and y-directionally varying periodic binary function 
%with the period of Tx and Ty, and  the width of wx and wy. 

%This function actually calculates the 2D Fourier series of the 2D periodic
%rectangle fuction.


Grating_gen_BinaryGrating_along_x;     % Fouirer series of grating structure along x
alpha_tm=0;
beta_tm=0;

% Grating_gen_BinaryGrating_along_y;     % Fourier series of grating structure along y
% alpha_tm=1;
% beta_tm=1;

grating_thick=Height/Nlay;               % grating thickness

lay_thick=zeros(Nlay,1); %  layer
for lnt=1:Nlay
lay_thick(lnt)=grating_thick; 
end;

ac_thick=zeros(1,Nlay); 
for cnt1=1:Nlay
	 cnt2=1;
    while cnt2 <= cnt1
       ac_thick(cnt1)=ac_thick(cnt1)+lay_thick(cnt2);
       cnt2=cnt2+1;
    end;
end; % for cnt1

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
% xx=5*[-Tx/2:Tx*0.01:Tx/2];
% yy=5*[-Ty/2:Ty*0.01:Ty/2];
% 
%      Gr_str=zeros(length(xx),length(yy));    
%      [xa ya]=meshgrid(xx,yy);
%      
%      for k=-2*nx:2*nx
%          for l=-2*ny:2*ny
%      
%              Gr_str=Gr_str+eps_xx(k+NBx,l+NBy,1)*exp(j*(k*xa*2*pi/Tx+l*ya*2*pi/Ty)); 
%          
%          end;
%      end;
%      
%      figure(5);mesh(yy,xx,abs(Gr_str)); xlabel('x-axis'); ylabel('y-axis');
%%
     
