%%
cell_num=7;

Nlay=NN*cell_num;    %  number of layers

epra=n_air^2  + i*1e-10;                       % air
eprd=n_diel^2 + i*1e-10;                       % dielectric

d=a;
dx=d/NN;
dz=dx;
sx=[-Tx/2+dx/2:dx:Tx/2-dx/2];   % x 방향 cell 중간좌표
% sz=cell_num*[-Tz/2+dz/2:dz:Tz/2-dz/2];   % z 방향 cell 중간좌표
sz=[-Tz/2+dz/2:dz:Tz/2-dz/2];   % z 방향 cell 중간좌표

Str=zeros(length(sx),length(sz));

% Layer A on the coordinate system (x,y,z)
Alay_thick=zeros(Nlay,1); %  layer
for k=1:Nlay
    Alay_thick(k)=dz;
end

Aac_thick=zeros(1,Nlay);
for cnt1=1:Nlay
    cnt2=1;
    while cnt2 <= cnt1
        Aac_thick(cnt1)=Aac_thick(cnt1)+Alay_thick(cnt2);
        cnt2=cnt2+1;
    end
end % for cnt1

% Grating_Gen_PCWG;
Grating_gen_PhC_resonator;
Aeps_L=eps_L;
Aaps_L=aps_L;
Amu_L=mu_L;
Abu_L=bu_L;

alpha_tm=1;
beta_tm=1;


%---------------------- Ae_h, Aa_h, Ag_h, Ab_h -----------------------
Aeps_xx=zeros(num_hx,num_hy,Nlay); Aeps_xy=zeros(num_hx,num_hy,Nlay); Aeps_xz=zeros(num_hx,num_hy,Nlay);   % permittivity
Aeps_yx=zeros(num_hx,num_hy,Nlay); Aeps_yy=zeros(num_hx,num_hy,Nlay); Aeps_yz=zeros(num_hx,num_hy,Nlay);
Aeps_zx=zeros(num_hx,num_hy,Nlay); Aeps_zy=zeros(num_hx,num_hy,Nlay); Aeps_zz=zeros(num_hx,num_hy,Nlay);

Aaps_xx=zeros(num_hx,num_hy,Nlay); Aaps_xy=zeros(num_hx,num_hy,Nlay); Aaps_xz=zeros(num_hx,num_hy,Nlay);   % reciprocal permittivity
Aaps_yx=zeros(num_hx,num_hy,Nlay); Aaps_yy=zeros(num_hx,num_hy,Nlay); Aaps_yz=zeros(num_hx,num_hy,Nlay);
Aaps_zx=zeros(num_hx,num_hy,Nlay); Aaps_zy=zeros(num_hx,num_hy,Nlay); Aaps_zz=zeros(num_hx,num_hy,Nlay);

Amu_xx=zeros(num_hx,num_hy,Nlay);  Amu_xy=zeros(num_hx,num_hy,Nlay);  Amu_xz=zeros(num_hx,num_hy,Nlay);    % permeability
Amu_yx=zeros(num_hx,num_hy,Nlay);  Amu_yy=zeros(num_hx,num_hy,Nlay);  Amu_yz=zeros(num_hx,num_hy,Nlay);
Amu_zx=zeros(num_hx,num_hy,Nlay);  Amu_zy=zeros(num_hx,num_hy,Nlay);  Amu_zz=zeros(num_hx,num_hy,Nlay);

Abu_xx=zeros(num_hx,num_hy,Nlay);  Abu_xy=zeros(num_hx,num_hy,Nlay);  Abu_xz=zeros(num_hx,num_hy,Nlay);    % reciprocal permeability
Abu_yx=zeros(num_hx,num_hy,Nlay);  Abu_yy=zeros(num_hx,num_hy,Nlay);  Abu_yz=zeros(num_hx,num_hy,Nlay);
Abu_zx=zeros(num_hx,num_hy,Nlay);  Abu_zy=zeros(num_hx,num_hy,Nlay);  Abu_zz=zeros(num_hx,num_hy,Nlay);

Aeps_xx=Aeps_L;
Aeps_yy=Aeps_L;
Aeps_zz=Aeps_L;
Aaps_xx=Aaps_L;
Aaps_yy=Aaps_L;
Aaps_zz=Aaps_L;
Amu_xx=  Amu_L;
Amu_yy=  Amu_L;
Amu_zz=  Amu_L;
Abu_xx=  Abu_L;
Abu_yy=  Abu_L;
Abu_zz=  Abu_L;

% figure(31);mesh(abs(Gr_str)); axis equal
%%
xx=1*(-Tx/2:Tx*0.002:Tx/2);
% yy=1*(-Ty/2:Ty*0.02:Ty/2);
yy=0;
zz=linspace(0,Tz,Nlay);

Gr_str=zeros(length(xx),Nlay);
[ya xa]=meshgrid(yy,xx);

for nn=1:Nlay
    Gr_str_=zeros(length(xx),length(yy));
    for k=-2*nx:2*nx
        for l=-2*ny:2*ny
            Gr_str_=Gr_str_+Aeps_xx(k+NBx,l+NBy,nn)*exp(j*(k*xa*2*pi/Tx+l*ya*2*pi/Ty));
        end
    end
    Gr_str(:,nn)=Gr_str_(:,1);
end

figure(62);set(gca,'fontsize',16);set(gca,'fontname','times new roman');
imagesc(zz/a,xx/a,sqrt(abs(repmat(Gr_str,[1 1]))));
set(gca,'ydir','normal');
set(gca,'fontname','times new roman');xlabel('z (in unit of a)');ylabel('x (in unit of a)');
axis([zz(1)/a zz(end)/a xx(1)/a xx(end)/a]);
set(gca,'fontname','times new roman');colorbar;caxis([0 4]);set(gca,'fontsize',16);


%%
% for visualization of Ez & Hz
L=NBx*NBy;
AEzz=zeros(L,L,Nlay);  %% permittivity
AAzz=zeros(L,L,Nlay);
AGzz=zeros(L,L,Nlay);
ABzz=zeros(L,L,Nlay);

teps_zz=zeros(num_hx,num_hy);
taps_zz=zeros(num_hx,num_hy);
tmu_zz=zeros(num_hx,num_hy);
tbu_zz=zeros(num_hx,num_hy);

tEzz=zeros(L,L);
tAzz=zeros(L,L);
tGzz=zeros(L,L);
tBzz=zeros(L,L);

for laynt=1:Nlay

    teps_zz=Aeps_zz(:,:,laynt);
    taps_zz=Aaps_zz(:,:,laynt);
    tmu_zz =Amu_zz(:,:,laynt);
    tbu_zz =Abu_zz(:,:,laynt);

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

    AEzz(:,:,laynt)=tEzz;
    AAzz(:,:,laynt)=tAzz;
    AGzz(:,:,laynt)=tGzz;
    ABzz(:,:,laynt)=tBzz;

end % for laynt



% Layer B  on the coordinate system (x',y',z')
%
%---------------------- Be_h, Ba_h, Bg_h, Bb_h -----------------------
Beps_xx=zeros(num_hx,num_hy,Nlay); Beps_xy=zeros(num_hx,num_hy,Nlay); Beps_xz=zeros(num_hx,num_hy,Nlay);   % permittivity
Beps_yx=zeros(num_hx,num_hy,Nlay); Beps_yy=zeros(num_hx,num_hy,Nlay); Beps_yz=zeros(num_hx,num_hy,Nlay);
Beps_zx=zeros(num_hx,num_hy,Nlay); Beps_zy=zeros(num_hx,num_hy,Nlay); Beps_zz=zeros(num_hx,num_hy,Nlay);

Baps_xx=zeros(num_hx,num_hy,Nlay); Baps_xy=zeros(num_hx,num_hy,Nlay); Baps_xz=zeros(num_hx,num_hy,Nlay);   % reciprocal permittivity
Baps_yx=zeros(num_hx,num_hy,Nlay); Baps_yy=zeros(num_hx,num_hy,Nlay); Baps_yz=zeros(num_hx,num_hy,Nlay);
Baps_zx=zeros(num_hx,num_hy,Nlay); Baps_zy=zeros(num_hx,num_hy,Nlay); Baps_zz=zeros(num_hx,num_hy,Nlay);

Bmu_xx=zeros(num_hx,num_hy,Nlay);  Bmu_xy=zeros(num_hx,num_hy,Nlay);  Bmu_xz=zeros(num_hx,num_hy,Nlay);    % permeability
Bmu_yx=zeros(num_hx,num_hy,Nlay);  Bmu_yy=zeros(num_hx,num_hy,Nlay);  Bmu_yz=zeros(num_hx,num_hy,Nlay);
Bmu_zx=zeros(num_hx,num_hy,Nlay);  Bmu_zy=zeros(num_hx,num_hy,Nlay);  Bmu_zz=zeros(num_hx,num_hy,Nlay);

Bbu_xx=zeros(num_hx,num_hy,Nlay);  Bbu_xy=zeros(num_hx,num_hy,Nlay);  Bbu_xz=zeros(num_hx,num_hy,Nlay);    % reciprocal permeability
Bbu_yx=zeros(num_hx,num_hy,Nlay);  Bbu_yy=zeros(num_hx,num_hy,Nlay);  Bbu_yz=zeros(num_hx,num_hy,Nlay);
Bbu_zx=zeros(num_hx,num_hy,Nlay);  Bbu_zy=zeros(num_hx,num_hy,Nlay);  Bbu_zz=zeros(num_hx,num_hy,Nlay);


Blay_thick=Alay_thick;
Bac_thick=Aac_thick;

Beps_L=Aeps_L;
Baps_L=Aaps_L;
Bmu_L=Amu_L;
Bbu_L=Abu_L;

Beps_xx=  Beps_L;
Beps_yy=  Beps_L;
Beps_zz=  Beps_L;
Baps_xx=  Baps_L;
Baps_yy=  Baps_L;
Baps_zz=  Baps_L;
Bmu_xx=  Bmu_L;
Bmu_yy=  Bmu_L;
Bmu_zz=  Bmu_L;
Bbu_xx=  Bbu_L;
Bbu_yy=  Bbu_L;
Bbu_zz=  Bbu_L;


L=NBx*NBy;
BEzz=zeros(L,L,Nlay);  %% permittivity
BAzz=zeros(L,L,Nlay);
BGzz=zeros(L,L,Nlay);
BBzz=zeros(L,L,Nlay);

teps_zz=zeros(num_hx,num_hy);
taps_zz=zeros(num_hx,num_hy);
tmu_zz=zeros(num_hx,num_hy);
tbu_zz=zeros(num_hx,num_hy);

tEzz=zeros(L,L);
tAzz=zeros(L,L);
tGzz=zeros(L,L);
tBzz=zeros(L,L);

for laynt=1:Nlay

    teps_zz=Beps_zz(:,:,laynt);
    taps_zz=Baps_zz(:,:,laynt);
    tmu_zz =Bmu_zz(:,:,laynt);
    tbu_zz =Bbu_zz(:,:,laynt);

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

    BEzz(:,:,laynt)=tEzz;
    BAzz(:,:,laynt)=tAzz;
    BGzz(:,:,laynt)=tGzz;
    BBzz(:,:,laynt)=tBzz;

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


