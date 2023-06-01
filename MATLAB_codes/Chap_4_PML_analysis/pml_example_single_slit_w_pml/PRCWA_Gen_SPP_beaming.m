% PRCWA_Gen_SPP_beaming.m
% Layer structure ------------------------------------------
Nlay=3;
% Nlay=1;

% material permittivity
epra=1;                       % air
eprm=(fun_Ag_nk(lambda))^2;   % silver

film_thick=300*nm;              % metal film thickness
% grating_thick=80*nm;            % dielectric grating thickness
grating_thick=0*nm;            % dielectric grating thickness

Grating_gen_SPP_beaming;        % Fourier series of grating structure
alpha_tm=1;
beta_tm=0;

lay_thick=zeros(Nlay,1);        %   layer
lay_thick(1)=film_thick;        %   metal substrate
lay_thick(2)=grating_thick;    %   dielectric grating
lay_thick(3)=0;                 %   half-infinite freespace with PML

ac_thick=zeros(1,Nlay); % ?????? ???? ?? ?????? ?????? ????
for cnt1=1:Nlay
    cnt2=1;
    while cnt2 <= cnt1
        ac_thick(cnt1)=ac_thick(cnt1)+lay_thick(cnt2);
        cnt2=cnt2+1;
    end
end % for cnt1

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

end % for lnt

%% Permittivity profile test
% xx=5*(-Tx/2:Tx*0.001:Tx/2)';
xx=1*(-Tx/2:Tx*0.001:Tx/2)';
% xx=0.07*(-Tx/2:Tx*0.001:Tx/2)';
Gr_str=zeros(length(xx),Nlay);

for laynt=1:Nlay
    for k=-2*nx:2*nx

        Gr_str(:,laynt)=Gr_str(:,laynt)+eps_xx(k+NBx,laynt)*exp(j*(k*xx*2*pi/Tx));

    end
end

figure(50);imagesc(abs(Gr_str));colorbar;

%% permittivity profile test
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
%          end
%      end
%
%      figure(5);mesh(abs(Gr_str));
%%


% for visualization of Ez & Hz
L=NBx*NBy;
Ezz=zeros(L,L,Nlay);  %% permittivity
Azz=zeros(L,L,Nlay);
Gzz=zeros(L,L,Nlay);
Bzz=zeros(L,L,Nlay);

teps_zz=zeros(num_hx,num_hy);
taps_zz=zeros(num_hx,num_hy);
tmu_zz=zeros(num_hx,num_hy);
tbu_zz=zeros(num_hx,num_hy);

tEzz=zeros(L,L);
tAzz=zeros(L,L);
tGzz=zeros(L,L);
tBzz=zeros(L,L);

for laynt=1:Nlay

    teps_zz=eps_zz(:,:,laynt);
    taps_zz=aps_zz(:,:,laynt);
    tmu_zz =mu_zz(:,:,laynt);
    tbu_zz =bu_zz(:,:,laynt);

    for k=1:NBx
        for l=1:NBy
            od_ind1=(k-1)*NBy+l;
            for kk=1:NBx
                for ll=1:NBy
                    od_ind2=(kk-1)*NBy+ll;

                    % permittivity
                    tEzz(od_ind1,od_ind2)=teps_zz(k-kk+NBx,l-ll+NBy);
                    tAzz(od_ind1,od_ind2)=taps_zz(k-kk+NBx,l-ll+NBy);
                    tGzz(od_ind1,od_ind2)=tmu_zz(k-kk+NBx,l-ll+NBy);
                    tBzz(od_ind1,od_ind2)=tbu_zz(k-kk+NBx,l-ll+NBy);
                end
            end
        end
    end

    Ezz(:,:,laynt)=tEzz;
    Azz(:,:,laynt)=tAzz;
    Gzz(:,:,laynt)=tGzz;
    Bzz(:,:,laynt)=tBzz;

end % for laynt


Kx=zeros(L,L);
Ky=zeros(L,L);
for k=1:NBx
    for l=1:NBy
        od_ind=(k-1)*NBy+l;
        Kx(od_ind,od_ind)=kx_vc(k)/k0;
        Ky(od_ind,od_ind)=ky_vc(l)/k0;
    end
end

