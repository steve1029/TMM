%----------------
sx=1+i*0.01;
sy=1;
sz=1;

pml_width=3*lambda;
pml_N=30;
gam=1/(1-i);

eps_surr=1;
muu_surr=1;

%----------------

% nx=100;NBx=2*nx+1;num_hx=4*nx+1;
% ny=0;NBy=2*ny+1;num_hy=4*ny+1;
% nx=10;NBx=2*nx+1;num_hx=4*nx+1;
% ny=0;NBy=2*ny+1;num_hy=4*ny+1;
% pml_N=3;

%% region set
d_x=Tx;
q_x=2*pml_width;
e_x=Tx-q_x;
% d_y=Ty;
% q_y=2*pml_width;
% e_y=Ty-q_y;

% --------------------------   +Ty/2
%      |             |
%      |             |
%   2  |      0      |  1
%      |             |
%      |             |
% --------------------------   -Ty/2
%   -e_x/2        +e_x/2  +d_x/2

Epsr_PML_xx=zeros(num_hx,num_hy);
Epsr_PML_yy=zeros(num_hx,num_hy);
Epsr_PML_zz=zeros(num_hx,num_hy);

Apsr_PML_xx=zeros(num_hx,num_hy);
Apsr_PML_yy=zeros(num_hx,num_hy);
Apsr_PML_zz=zeros(num_hx,num_hy);

Mpsr_PML_xx=zeros(num_hx,num_hy);
Mpsr_PML_yy=zeros(num_hx,num_hy);
Mpsr_PML_zz=zeros(num_hx,num_hy);

Bpsr_PML_xx=zeros(num_hx,num_hy);
Bpsr_PML_yy=zeros(num_hx,num_hy);
Bpsr_PML_zz=zeros(num_hx,num_hy);


rect_region0=zeros(num_hx,num_hy);% region 0

%% region 1-2

pml_num=(1:pml_N);
pml_x=(Tx/2-pml_width) + pml_width*([pml_num pml_N+1]-1)/pml_N;
% pml_y=(Ty/2-pml_width) + pml_width*([pml_num pml_N+1]-1)/paml_N;


for mm=1:num_hx
%     disp([num2str(mm) ' / ' num2str(num_hx)]);
    m=mm-NBx;
    for nn=1:num_hy
        n=nn-NBy;
        
        % region 0
        rect_region0(mm,nn)=rect_2D(m,n,1,Tx,Ty,-e_x/2,e_x/2,-Ty/2,Ty/2);

        for k=1:pml_N-0  % 

            % region 1 and 2
            pml_x1=pml_x(k);
            pml_x2=pml_x(k+1);
            cos2sx=(1-gam)*cos(pi * (abs(pml_x1)-e_x/2) / q_x)^2/sx;
            cos2sy=1;
            rect_form=rect_2D(m,n,1,Tx,Ty, pml_x1, pml_x2,-Ty/2,Ty/2);
            PML_stretched_evaluate;
            rect_form=rect_2D(m,n,1,Tx,Ty,-pml_x2,-pml_x1,-Ty/2,Ty/2);
            PML_stretched_evaluate;

        end
    end
end



%%
xx=1*(-Tx/2:Tx*0.01:Tx/2);
yy=0;
% yy=1*(-Ty/2:Ty*0.01:Ty/2);

Gr_strEx=zeros(length(xx),length(yy));
Gr_strEy=zeros(length(xx),length(yy));
Gr_strEz=zeros(length(xx),length(yy));
Gr_strAx=zeros(length(xx),length(yy));
Gr_strAy=zeros(length(xx),length(yy));
Gr_strAz=zeros(length(xx),length(yy));
[ya xa]=meshgrid(yy,xx);

for k=-2*nx:2*nx
    for l=-2*ny:2*ny
        aa=exp(j*(k*xa*2*pi/Tx+l*ya*2*pi/Ty));
        Gr_strEx=Gr_strEx+Epsr_PML_xx(k+NBx,l+NBy)*aa;
        Gr_strEy=Gr_strEy+Epsr_PML_yy(k+NBx,l+NBy)*aa;
        Gr_strEz=Gr_strEz+Epsr_PML_zz(k+NBx,l+NBy)*aa;
        Gr_strAx=Gr_strAx+Apsr_PML_xx(k+NBx,l+NBy)*aa;
        Gr_strAy=Gr_strAy+Apsr_PML_yy(k+NBx,l+NBy)*aa;
        Gr_strAz=Gr_strAz+Apsr_PML_zz(k+NBx,l+NBy)*aa;
    end
end

figure(10);
subplot(2,3,1);imagesc(yy/um,xx/um,real(Gr_strEx));title('epsx');set(gca,'ydir','normal');colorbar;
subplot(2,3,2);imagesc(yy/um,xx/um,real(Gr_strEy));title('epsy');set(gca,'ydir','normal');colorbar;
subplot(2,3,3);imagesc(yy/um,xx/um,real(Gr_strEz));title('epsz');set(gca,'ydir','normal');colorbar;
subplot(2,3,4);imagesc(yy/um,xx/um,real(Gr_strAx));title('apsx');set(gca,'ydir','normal');colorbar;
subplot(2,3,5);imagesc(yy/um,xx/um,real(Gr_strAy));title('apsy');set(gca,'ydir','normal');colorbar;
subplot(2,3,6);imagesc(yy/um,xx/um,real(Gr_strAz));title('apsz');set(gca,'ydir','normal');colorbar;





