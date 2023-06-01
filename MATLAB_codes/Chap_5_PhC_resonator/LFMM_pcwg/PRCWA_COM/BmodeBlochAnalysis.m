% B mode Bloch Mode Analysis

L=NBx*NBy;

Ta=zeros(2*L,2*L,Nlay); % left to rignt
Ra=zeros(2*L,2*L,Nlay); % left to right
Tb=zeros(2*L,2*L,Nlay); % right to left
Rb=zeros(2*L,2*L,Nlay); % right to left
Ca=zeros(4*L,2*L,Nlay); % left to right
Cb=zeros(4*L,2*L,Nlay); % left to right
tCa=zeros(4*L,2*L,Nlay); % left to right
tCb=zeros(4*L,2*L,Nlay); % left to right

%% Layer S-matrix computation of layers

eps_xx=Beps_xx; eps_xy=Beps_xy; eps_xz=Beps_xz;
eps_yx=Beps_yx; eps_yy=Beps_yy; eps_yz=Beps_yz;
eps_zx=Beps_zx; eps_zy=Beps_zy; eps_zz=Beps_zz;

aps_xx=Baps_xx; aps_xy=Baps_xy; aps_xz=Baps_xz;
aps_yx=Baps_yx; aps_yy=Baps_yy; aps_yz=Baps_yz;
aps_zx=Baps_zx; aps_zy=Baps_zy; aps_zz=Baps_zz;

mu_xx=Bmu_xx; mu_xy=Bmu_xy; mu_xz=Bmu_xz;
mu_yx=Bmu_yx; mu_yy=Bmu_yy; mu_yz=Bmu_yz;
mu_zx=Bmu_zx; mu_zy=Bmu_zy; mu_zz=Bmu_zz;

bu_xx=Bbu_xx; bu_xy=Bbu_xy; bu_xz=Bbu_xz;
bu_yx=Bbu_yx; bu_yy=Bbu_yy; bu_yz=Bbu_yz;
bu_zx=Bbu_zx; bu_zy=Bbu_zy; bu_zz=Bbu_zz;

lay_thick=Blay_thick;

%Off_diagonal_tensor_SMM;   
Diagonal_tensor_SMM;

%% Layer S-matrix of grating body


   I=eye(2*L,2*L);
   T_temp1a=Ta(:,:,1);
   R_temp1a=Ra(:,:,1);
   T_temp1b=Tb(:,:,1);
   R_temp1b=Rb(:,:,1);
   
   %%% Important
   tCa=Ca;
   tCb=Cb;
   %%%
   
for laynt=2:Nlay
   
   %laynt
   
   T_temp2a=Ta(:,:,laynt);
   R_temp2a=Ra(:,:,laynt);
   T_temp2b=Tb(:,:,laynt);
   R_temp2b=Rb(:,:,laynt);

	RRa=(R_temp1a+T_temp1b*inv(I-R_temp2a*R_temp1b)*R_temp2a*T_temp1a);
	TTa=T_temp2a*inv(I-R_temp1b*R_temp2a)*T_temp1a;

	RRb=(R_temp2b+T_temp2a*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b);
   TTb=T_temp1b*inv(I-R_temp2a*R_temp1b)*T_temp2b;
   
   
   for k=1:laynt-1
   	
   tCa(:,:,k)=Ca(:,:,k)+Cb(:,:,k)*inv(I-R_temp2a*R_temp1b)*R_temp2a*T_temp1a;
   tCb(:,:,k)=Cb(:,:,k)*inv(I-R_temp2a*R_temp1b)*T_temp2b;

   end; % for k
   
   tCa(:,:,laynt)=Ca(:,:,laynt)*inv(I-R_temp1b*R_temp2a)*T_temp1a;
   tCb(:,:,laynt)=Cb(:,:,laynt)+Ca(:,:,laynt)*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b;
   
   
   
   T_temp1a=TTa;
   R_temp1a=RRa;
   T_temp1b=TTb;
   R_temp1b=RRb;
   
   Ca=tCa;
   Cb=tCb; 
end; % laynt

   TTa=T_temp1a;  % left-to-right transmission operator     
   RRa=R_temp1a;  % left-to-right reflection operator
   TTb=T_temp1b;  % right-to-left transmission operator
   RRb=R_temp1b;  % right-to-left reflection operator

  
%%   positive mode
SBSa=[TTa zeros(2*L,2*L); RRa -I];
SBSb=[I   -RRb;     zeros(2*L,2*L)  -TTb];

[W,D]=eig(SBSa,SBSb);

for k=1:4*L
e_value(k)=log(D(k,k))/d;
end;

pnl=0;
mnl=0;
for k=1:4*L
    if abs(real(e_value(k))) >= abs(imag(e_value(k)))
        pnl=pnl+1;pvalue(pnl)=e_value(k);pvec(:,pnl)=W(:,k);
    else
        mnl=mnl+1;mvalue(mnl)=e_value(k);mvec(:,mnl)=W(:,k);
    end;
end;

if pnl ~=0 
[pevv pevn]=sort(abs(real(pvalue)));

for k=1:pnl
e_value(k)=pvalue(pevn(k));
e_vector(:,k)=pvec(:,pevn(k));
end;
end;

[mevv mevn]=sort(abs(imag(mvalue)));
for k=1:mnl
    e_value(k+pnl)=mvalue(mevn(k));
    e_vector(:,k+pnl)=mvec(:,mevn(k));

end;

pcnt=0;
%mcnt=0;

for k=1:4*L
   if real(e_value(k)) <=0
      pcnt=pcnt+1;  % positive mode
   else 
 %  	mcnt=mcnt+1;    % negative mode   
   end;   
end;


pe_value=zeros(1,pcnt); 			% positive mode  real(e_value)<0
%me_value=zeros(1,mcnt);				% negative mode  real(e_value)>0
pe_vector=zeros(4*L,pcnt);			% positive mode  eigen_vector
%me_vector=zeros(4*L,mcnt);			% negative mode  eigen_vector

pcnt=0;
%mcnt=0;

for k=1:4*L
   
	if real(e_value(k)) <= 0   
      pcnt=pcnt+1; 
      pe_value(pcnt)=e_value(k);
      pe_vector(:,pcnt)=e_vector(:,k);
   else
 %     mcnt=mcnt+1;
 %     me_value(mcnt)=e_value(k);
 %     me_vector(:,mcnt)=e_vector(:,k);
   end;
   
end;


%%  negative mode

[W, D]=eig(SBSb,SBSa);

for k=1:4*L
e_value(k)=-log(D(k,k))/d;
end;

pnl=0;
mnl=0;
for k=1:4*L
    if abs(real(e_value(k))) >= abs(imag(e_value(k)))
        pnl=pnl+1;pvalue(pnl)=e_value(k);pvec(:,pnl)=W(:,k);
    else
      mnl=mnl+1;mvalue(mnl)=e_value(k);mvec(:,mnl)=W(:,k);
    end;
end;


if pnl ~= 0
[pevv pevn]=sort(abs(real(pvalue)));

for k=1:pnl
e_value(k)=pvalue(pevn(k));
e_vector(:,k)=pvec(:,pevn(k));
end;
end;
[mevv mevn]=sort(abs(imag(mvalue)));
for k=1:mnl
    e_value(k+pnl)=mvalue(mevn(k));
e_vector(:,k+pnl)=mvec(:,mevn(k));

end;

%pcnt=0;
mcnt=0;

for k=1:4*L
   if real(e_value(k)) <=0
%      pcnt=pcnt+1;  % positive mode
   else 
   	mcnt=mcnt+1;    % negative mode   
   end;   
end;


%pe_value=zeros(1,pcnt); 			% positive mode  real(e_value)<0
me_value=zeros(1,mcnt);				% negative mode  real(e_value)>0
%pe_vector=zeros(4*L,pcnt);			% positive mode  eigen_vector
me_vector=zeros(4*L,mcnt);			% negative mode  eigen_vector

%pcnt=0;
mcnt=0;

for k=1:4*L
   
	if real(e_value(k)) <= 0   
 %     pcnt=pcnt+1; 
 %     pe_value(pcnt)=e_value(k);
 %     pe_vector(:,pcnt)=e_vector(:,k);
   else
      mcnt=mcnt+1;
      me_value(mcnt)=e_value(k);
      me_vector(:,mcnt)=e_vector(:,k);
   end;
   
end;




CCp=zeros(4*L,Nlay,pcnt); % positive Bloch mode
CCn=zeros(4*L,Nlay,mcnt); % negative Bloch mode
 
        for k=1:pcnt
           for lnt=1:Nlay
               CCp(:,lnt,k)=Ca(:,:,lnt)*pe_vector(1:2*L,k)+Cb(:,:,lnt)*exp(pe_value(k)*d)*pe_vector(2*L+1:4*L,k);
            end;
        end;
        
       for k=1:mcnt
           for lnt=1:Nlay
               CCn(:,lnt,k)=Ca(:,:,lnt)*me_vector(1:2*L,k)+Cb(:,:,lnt)*exp(me_value(k)*d)*me_vector(2*L+1:4*L,k);
            end;
        end;
               
%%  pseudo-Fourier representation of Bloch modes        
%         
% nz=nx;
% NBz=2*nz+1;
% num_hz=2*NBz-1; 

Tz=d;
dz=Tz/NBz;
zz=[0:dz:Tz-dz];
zn=length(zz);

rg2cnt=zeros(Nlay,1);   % region 2 : grating region 

for znt=1:zn
    
    zv=zz(znt);
    
    if (0 <= zv) & (zv <= Bac_thick(Nlay))
            
       for laynt=1:Nlay
             if  ( ( Bac_thick(laynt)-Blay_thick(laynt) ) < zv ) & ( zv <= Bac_thick(laynt)  )
              rg2cnt(laynt)=rg2cnt(laynt)+1; 
             end;
       end;
    
          if zv==0
            rg2cnt(1)=rg2cnt(1)+1; 
       end;
       
    end;
       
    
end; % for znt


% positive Bloch Mode ----------------------------------------------
%
 pEFy=zeros(NBx,NBy,NBz,pcnt);     
 pEFx=zeros(NBx,NBy,NBz,pcnt);
 pEFz=zeros(NBx,NBy,NBz,pcnt);
 
 pHFy=zeros(NBx,NBy,NBz,pcnt);     
 pHFx=zeros(NBx,NBy,NBz,pcnt);
 pHFz=zeros(NBx,NBy,NBz,pcnt);
 
 % region 2 : grating region
 for laynt=1:Nlay
     if laynt==1
 zz2ind(laynt,1)=1;
 zz2ind(laynt,2)=zz2ind(laynt,1)+rg2cnt(laynt)-1 ;
 
     else
 zz2ind(laynt,1)=zz2ind(laynt-1,2)+1;
 zz2ind(laynt,2)=zz2ind(laynt,1)+rg2cnt(laynt)-1 ;
     end;
 end;
 
 fc_ind=zz2ind;  % field calculation index for region 2

 
for pbm_ind=1:pcnt
 
  C_p=CCp(1:2*L,    :,pbm_ind);
  C_n=CCp(2*L+1:4*L,:,pbm_ind);
  
 Ex_z=zeros(NBx,NBy,NBz);
 Ey_z=zeros(NBx,NBy,NBz);     
 Ez_z=zeros(NBx,NBy,NBz);     
 Hx_z=zeros(NBx,NBy,NBz);
 Hy_z=zeros(NBx,NBy,NBz);     
 Hz_z=zeros(NBx,NBy,NBz);
 
  Iden=zeros(2*L,2*L);
  for k=1:2*L
    Iden(k,k)=1;   
  end;
 
   
 for laynt=1:Nlay
    
    pevalue=Peigvalue(1,:,laynt);
    mevalue=Meigvalue(1,:,laynt);

    pcnt=length(pevalue); % number of positive modes
    mcnt=length(mevalue); % number of negative modes

    zm=Bac_thick(laynt)-Blay_thick(laynt);
    zp=Bac_thick(laynt);

    zzp=zz(zz2ind(laynt,1):zz2ind(laynt,2))-zm;
    zzm=zz(zz2ind(laynt,1):zz2ind(laynt,2))-zp;
    
  for zz_cnt=1:length(zzp)
   
    for k=1:2*L
        X(k,k)   =exp(pevalue(k)*zzp(zz_cnt));
        Xinv(k,k)=exp(mevalue(k)*zzm(zz_cnt));
    end;
    
    tmp1=X*C_p(:,laynt)*exp(-pe_value(pbm_ind)*(zzp(zz_cnt)+zm));
    pfEx=Pf_Ex(:,:,laynt)*tmp1;
    pfEy=Pf_Ey(:,:,laynt)*tmp1;
    pfEz=Pf_Ez(:,:,laynt)*tmp1;
    pfHx=Pf_Hx(:,:,laynt)*tmp1;
    pfHy=Pf_Hy(:,:,laynt)*tmp1;
    pfHz=Pf_Hz(:,:,laynt)*tmp1;
    
    tmp2=Xinv*C_n(:,laynt)*exp(-pe_value(pbm_ind)*(zzp(zz_cnt)+zm));
    mfEx=Mf_Ex(:,:,laynt)*tmp2;
    mfEy=Mf_Ey(:,:,laynt)*tmp2;
    mfEz=Mf_Ez(:,:,laynt)*tmp2;
    mfHx=Mf_Hx(:,:,laynt)*tmp2;
    mfHy=Mf_Hy(:,:,laynt)*tmp2;
    mfHz=Mf_Hz(:,:,laynt)*tmp2;
    
    
     for k=1:NBx
      for l=1:NBy
          
          Ex_z(k,l,zz_cnt+fc_ind(laynt,1)-1)= pfEx((k-1)*NBy+l) + mfEx((k-1)*NBy+l) ;   % Ex_xz    
          Ey_z(k,l,zz_cnt+fc_ind(laynt,1)-1)= pfEy((k-1)*NBy+l) + mfEy((k-1)*NBy+l) ;   % Ey_xz
          Ez_z(k,l,zz_cnt+fc_ind(laynt,1)-1)= pfEz((k-1)*NBy+l) + mfEz((k-1)*NBy+l) ;   % Ex_xz    
          
          Hx_z(k,l,zz_cnt+fc_ind(laynt,1)-1)= pfHx((k-1)*NBy+l) + mfHx((k-1)*NBy+l) ;   % Hx_xz
          Hy_z(k,l,zz_cnt+fc_ind(laynt,1)-1)= pfHy((k-1)*NBy+l) + mfHy((k-1)*NBy+l) ;   % Hy_xz
          Hz_z(k,l,zz_cnt+fc_ind(laynt,1)-1)= pfHz((k-1)*NBy+l) + mfHz((k-1)*NBy+l) ;   % Hx_xz

      end; % for l
  end; % for k
      
  end; % zz_cnt

end; % for laynt
 

%%% z-direction Fourier coefficients
ey=zeros(1,NBz);
ex=zeros(1,NBz);
ez=zeros(1,NBz);

hy=zeros(1,NBz);
hx=zeros(1,NBz);
hz=zeros(1,NBz);


for k=1:NBx
    for l=1:NBy
        
        for znt=1:NBz
        ey(znt)=Ey_z(k,l,znt);
        ex(znt)=Ex_z(k,l,znt);
        ez(znt)=Ez_z(k,l,znt);
        hy(znt)=Hy_z(k,l,znt);
        hx(znt)=Hx_z(k,l,znt);
        hz(znt)=Hz_z(k,l,znt);
        end;
        
        
        bb=DFT_conv(ex,NBz);
        for znt=1:NBz
        pEFx(k,l,znt,pbm_ind)=bb(znt);
        end;
        
        
        bb=DFT_conv(ey,NBz);
        for znt=1:NBz
        pEFy(k,l,znt,pbm_ind)=bb(znt);
        end;
        
        
        bb=DFT_conv(ez,NBz);
        for znt=1:NBz
        pEFz(k,l,znt,pbm_ind)=bb(znt);
        end;
        
        
        bb=DFT_conv(hx,NBz);
        for znt=1:NBz
        pHFx(k,l,znt,pbm_ind)=bb(znt);
        end;
        
        
        bb=DFT_conv(hy,NBz);
        for znt=1:NBz
        pHFy(k,l,znt,pbm_ind)=bb(znt);
        end;
        
        
        bb=DFT_conv(hz,NBz);
        for znt=1:NBz
        pHFz(k,l,znt,pbm_ind)=bb(znt);
        end;
        
        
    end;    % for l
end;    % for k

end; % for pbm_ind= 1 : pcnt  
        
% negative Bloch Mode -----------------------------------

 nEFy=zeros(NBx,NBy,NBz,mcnt);     
 nEFx=zeros(NBx,NBy,NBz,mcnt);
 nEFz=zeros(NBx,NBy,NBz,mcnt);
 
 nHFy=zeros(NBx,NBy,NBz,mcnt);     
 nHFx=zeros(NBx,NBy,NBz,mcnt);
 nHFz=zeros(NBx,NBy,NBz,mcnt);
 
 
for pbm_ind=1:mcnt

  C_p=CCn(1:2*L,    :,pbm_ind);
  C_n=CCn(2*L+1:4*L,:,pbm_ind);
  
  Ex_z=zeros(NBx,NBy,NBz);
  Ey_z=zeros(NBx,NBy,NBz);     
  Ez_z=zeros(NBx,NBy,NBz);     
  Hx_z=zeros(NBx,NBy,NBz);
  Hy_z=zeros(NBx,NBy,NBz);     
  Hz_z=zeros(NBx,NBy,NBz);

  Iden=zeros(2*L,2*L);
  for k=1:2*L
    Iden(k,k)=1;   
  end;
 
   
 for laynt=1:Nlay
    
    pevalue=Peigvalue(1,:,laynt);
    mevalue=Meigvalue(1,:,laynt);

    pcnt=length(pevalue); % number of positive modes
    mcnt=length(mevalue); % number of negative modes

    zm=Aac_thick(laynt)-Alay_thick(laynt);
    zp=Aac_thick(laynt);

    zzp=zz(zz2ind(laynt,1):zz2ind(laynt,2))-zm;
    zzm=zz(zz2ind(laynt,1):zz2ind(laynt,2))-zp;
    
  for zz_cnt=1:length(zzp)
   
    for k=1:2*L
        X(k,k)   =exp(pevalue(k)*zzp(zz_cnt));
        Xinv(k,k)=exp(mevalue(k)*zzm(zz_cnt));
    end;
    
    tmp1=X*C_p(:,laynt)*exp(-me_value(pbm_ind)*(zzm(zz_cnt)+zp));
    pfEx=Pf_Ex(:,:,laynt)*tmp1;
    pfEy=Pf_Ey(:,:,laynt)*tmp1;
    pfEz=Pf_Ez(:,:,laynt)*tmp1;
    pfHx=Pf_Hx(:,:,laynt)*tmp1;
    pfHy=Pf_Hy(:,:,laynt)*tmp1;
    pfHz=Pf_Hz(:,:,laynt)*tmp1;
    
    tmp2=Xinv*C_n(:,laynt)*exp(-me_value(pbm_ind)*(zzm(zz_cnt)+zp));
    mfEx=Mf_Ex(:,:,laynt)*tmp2;
    mfEy=Mf_Ey(:,:,laynt)*tmp2;
    mfEz=Mf_Ez(:,:,laynt)*tmp2;
    mfHx=Mf_Hx(:,:,laynt)*tmp2;
    mfHy=Mf_Hy(:,:,laynt)*tmp2;
    mfHz=Mf_Hz(:,:,laynt)*tmp2;
    
    
    
     for k=1:NBx
      for l=1:NBy
          
          Ex_z(k,l,zz_cnt+fc_ind(laynt,1)-1)= pfEx((k-1)*NBy+l) + mfEx((k-1)*NBy+l) ;   % Ex_xz    
          Ey_z(k,l,zz_cnt+fc_ind(laynt,1)-1)= pfEy((k-1)*NBy+l) + mfEy((k-1)*NBy+l) ;   % Ey_xz
          Ez_z(k,l,zz_cnt+fc_ind(laynt,1)-1)= pfEz((k-1)*NBy+l) + mfEz((k-1)*NBy+l) ;   % Ex_xz    
          
          Hx_z(k,l,zz_cnt+fc_ind(laynt,1)-1)= pfHx((k-1)*NBy+l) + mfHx((k-1)*NBy+l) ;   % Hx_xz
          Hy_z(k,l,zz_cnt+fc_ind(laynt,1)-1)= pfHy((k-1)*NBy+l) + mfHy((k-1)*NBy+l) ;   % Hy_xz
          Hz_z(k,l,zz_cnt+fc_ind(laynt,1)-1)= pfHz((k-1)*NBy+l) + mfHz((k-1)*NBy+l) ;   % Hx_xz

      end; % for l
  end; % for k
  
  end; % zz_cnt

end; % for laynt


%%% z-direction Fourier coefficients
ey=zeros(1,NBz);
ex=zeros(1,NBz);
ez=zeros(1,NBz);

hy=zeros(1,NBz);
hx=zeros(1,NBz);
hz=zeros(1,NBz);


for k=1:NBx
    for l=1:NBy
        
        for znt=1:NBz
        ey(znt)=Ey_z(k,l,znt);
        ex(znt)=Ex_z(k,l,znt);
        ez(znt)=Ez_z(k,l,znt);
        hy(znt)=Hy_z(k,l,znt);
        hx(znt)=Hx_z(k,l,znt);
        hz(znt)=Hz_z(k,l,znt);
        end;
        
        
        bb=DFT_conv(ex,NBz);
        for znt=1:NBz
        nEFx(k,l,znt,pbm_ind)=bb(znt);
        end;
        
        
        bb=DFT_conv(ey,NBz);
        for znt=1:NBz
        nEFy(k,l,znt,pbm_ind)=bb(znt);
        end;
        
        
        bb=DFT_conv(ez,NBz);
        for znt=1:NBz
        nEFz(k,l,znt,pbm_ind)=bb(znt);
        end;
        
        
        bb=DFT_conv(hx,NBz);
        for znt=1:NBz
        nHFx(k,l,znt,pbm_ind)=bb(znt);
        end;
        
        
        bb=DFT_conv(hy,NBz);
        for znt=1:NBz
        nHFy(k,l,znt,pbm_ind)=bb(znt);
        end;
        
        
        bb=DFT_conv(hz,NBz);
        for znt=1:NBz
        nHFz(k,l,znt,pbm_ind)=bb(znt);
        end;
        
        
    end;    % for l
end;    % for k

end; % for pbm_ind= 1 : mcnt  

% obtained result-------------------------------------------------------
%  
%  pe_value=zeros(1,pcnt); 			
%  me_value=zeros(1,mcnt);			
%  pEFy=zeros(NBx,NBy,NBz,pcnt);     
%  pEFx=zeros(NBx,NBy,NBz,pcnt);
%  pEFz=zeros(NBx,NBy,NBz,pcnt);
%  
%  pHFy=zeros(NBx,NBy,NBz,pcnt);     
%  pHFx=zeros(NBx,NBy,NBz,pcnt);
%  pHFz=zeros(NBx,NBy,NBz,pcnt);
% 
%  nEFy=zeros(NBx,NBy,NBz,mcnt);     
%  nEFx=zeros(NBx,NBy,NBz,mcnt);
%  nEFz=zeros(NBx,NBy,NBz,mcnt);
%  
%  nHFy=zeros(NBx,NBy,NBz,mcnt);     
%  nHFx=zeros(NBx,NBy,NBz,mcnt);
%  nHFz=zeros(NBx,NBy,NBz,mcnt);
%-------------------------------------------------------------------------