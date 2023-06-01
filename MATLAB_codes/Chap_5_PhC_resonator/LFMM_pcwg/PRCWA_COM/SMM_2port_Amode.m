% SMM_2port_Amode

Iden=zeros(2*L,2*L);
for k=1:2*L
Iden(k,k)=1;   
end;

KII=zeros(2*L,2*L); % KII matrix
for k=1:NBx
   for l=1:NBy
od_ind1=(k-1)*NBy+l;      
K11(od_ind1,od_ind1)=(Akx_vc(k)*Aky_vc(l))/Akz_vc(k,l);
K12(od_ind1,od_ind1)=(Akz_vc(k,l)^2+Akx_vc(k)^2)/Akz_vc(k,l);
K21(od_ind1,od_ind1)=-(Aky_vc(l)^2+Akz_vc(k,l)^2)/Akz_vc(k,l);
K22(od_ind1,od_ind1)=-( Aky_vc(l)*Akx_vc(k) )/Akz_vc(k,l);
   end;
end;

KII(1:L, 1:L)=K11;
KII(1:L, L+1:2*L)=K12;
KII(L+1:2*L, 1:L)=K21;
KII(L+1:2*L,L+1:2*L )=K22;

AWh=Iden;
AVh=KII*(1/(w0*mu0));

AWpm=zeros(2*L,pcnt); AVpm=zeros(2*L,pcnt); 
AWnm=zeros(2*L,mcnt); AVnm=zeros(2*L,mcnt);

AWpp=zeros(2*L,pcnt); AVpp=zeros(2*L,pcnt);
AWnp=zeros(2*L,mcnt); AVnp=zeros(2*L,mcnt);
            
for pbm_ind=1:pcnt % AWpm, AVpm, AWpp, AVpp

    for k=1:NBx
        for l=1:NBy
   
    [ey ex]=awp(zm-zm,k,l,pbm_ind,ApEFy,ApEFx,ap_evalue);
    AWpm(NBy*(k-1)+l,pbm_ind)=ey;
    AWpm(L+NBy*(k-1)+l,pbm_ind)=ex;
    
    [hy hx]=avp(zm-zm,k,l,pbm_ind,ApHFy,ApHFx,ap_evalue);
    AVpm(NBy*(k-1)+l,pbm_ind)=hy;
    AVpm(L+NBy*(k-1)+l,pbm_ind)=hx;

    
    [ey ex]=awp(zp-zm,k,l,pbm_ind,ApEFy,ApEFx,ap_evalue);
    AWpp(NBy*(k-1)+l,pbm_ind)=ey;
    AWpp(L+NBy*(k-1)+l,pbm_ind)=ex;
    
    [hy hx]=avp(zp-zm,k,l,pbm_ind,ApHFy,ApHFx,ap_evalue);
    AVpp(NBy*(k-1)+l,pbm_ind)=hy;
    AVpp(L+NBy*(k-1)+l,pbm_ind)=hx;
    
        end;
    end;

end;

for pbm_ind=1:mcnt % AWnm, AVnm, AWnp, AVnp
    
    for k=1:NBx
        for l=1:NBy
   
    [ey ex]=awn(zm-zp,k,l,pbm_ind,AnEFy,AnEFx,am_evalue);
    AWnm(NBy*(k-1)+l,pbm_ind)=ey;
    AWnm(L+NBy*(k-1)+l,pbm_ind)=ex;
    
    [hy hx]=avn(zm-zp,k,l,pbm_ind,AnHFy,AnHFx,am_evalue);
    AVnm(NBy*(k-1)+l,pbm_ind)=hy;
    AVnm(L+NBy*(k-1)+l,pbm_ind)=hx;

    
    [ey ex]=awn(zp-zp,k,l,pbm_ind,AnEFy,AnEFx,am_evalue);
    AWnp(NBy*(k-1)+l,pbm_ind)=ey;
    AWnp(L+NBy*(k-1)+l,pbm_ind)=ex;
    
    [hy hx]=avn(zp-zp,k,l,pbm_ind,AnHFy,AnHFx,am_evalue);
    AVnp(NBy*(k-1)+l,pbm_ind)=hy;
    AVnp(L+NBy*(k-1)+l,pbm_ind)=hx;
    
        end;
    end;

end;

% Layer S-matrix
U=eye(2*NBx*NBy);
S11=inv(AWh)*AWpm+inv(AVh)*AVpm;
S12=inv(AWh)*AWnm+inv(AVh)*AVnm;
S21=inv(AWh)*AWpp-inv(AVh)*AVpp;
S22=inv(AWh)*AWnp-inv(AVh)*AVnp;
S=[S11 S12;S21 S22];

% left-to-right layer S-matrix
ACa=inv(S)*[2*U;zeros(2*NBx*NBy)];   %% very good!!!
ACap=ACa(1:pcnt,:);
ACam=ACa(pcnt+1:pcnt+mcnt,:);
AR11=inv(AWh)*(AWpm*ACap+AWnm*ACam-AWh*U);
AT12=inv(AWh)*(AWpp*ACap+AWnp*ACam);

% right-to-left layer S-matrix
ACb=inv(S)*[zeros(2*NBx*NBy);2*U];
ACbp=ACb(1:pcnt,:);
ACbm=ACb(pcnt+1:pcnt+mcnt,:);
AR22=inv(AWh)*(AWpp*ACbp+AWnp*ACbm-AWh*U);
AT21=inv(AWh)*(AWpm*ACbp+AWnm*ACbm);

% Boundary S-matrix

AWp=zeros(2*L,pcnt); AVp=zeros(2*L,pcnt); 
AWn=zeros(2*L,mcnt); AVn=zeros(2*L,mcnt);


for pbm_ind=1:pcnt % AWpm, AVpm, AWpp, AVpp

    for k=1:NBx
        for l=1:NBy
   
    [ey ex]=awp(zc,k,l,pbm_ind,ApEFy,ApEFx,ap_evalue);
    AWp(NBy*(k-1)+l,pbm_ind)=ey;
    AWp(L+NBy*(k-1)+l,pbm_ind)=ex;
    
    [hy hx]=avp(zc,k,l,pbm_ind,ApHFy,ApHFx,ap_evalue);
    AVp(NBy*(k-1)+l,pbm_ind)=hy;
    AVp(L+NBy*(k-1)+l,pbm_ind)=hx;

    
        end;
    end;

end;

for pbm_ind=1:mcnt % AWnm, AVnm, AWnp, AVnp
    
    for k=1:NBx
        for l=1:NBy
   
    [ey ex]=awn(zc,k,l,pbm_ind,AnEFy,AnEFx,am_evalue);
    AWn(NBy*(k-1)+l,pbm_ind)=ey;
    AWn(L+NBy*(k-1)+l,pbm_ind)=ex;
    
    [hy hx]=avn(zc,k,l,pbm_ind,AnHFy,AnHFx,am_evalue);
    AVn(NBy*(k-1)+l,pbm_ind)=hy;
    AVn(L+NBy*(k-1)+l,pbm_ind)=hx;

    
        end;
    end;

end;

% left side semi-infinite waveguide
ALR11=-inv(inv(AWh)*AWn-inv(AVh)*AVn)*(inv(AWh)*AWp-inv(AVh)*AVp);
ALT12=inv(inv(AWn)*AWh-inv(AVn)*AVh)*(inv(AWn)*AWp-inv(AVn)*AVp);
ALT21=2*inv(inv(AWh)*AWn-inv(AVh)*AVn);
ALR22=-inv(inv(AWn)*AWh-inv(AVn)*AVh)*(inv(AWn)*AWh+inv(AVn)*AVh);

% right side semi-infinite waveguide
ARR11=-inv(inv(AWp)*AWh+inv(AVp)*AVh)*(inv(AWp)*AWh-inv(AVp)*AVh);
ART12=2*inv(inv(AWh)*AWp+inv(AVh)*AVp);
ART21=inv(inv(AWp)*AWh+inv(AVp)*AVh)*(inv(AWp)*AWn-inv(AVp)*AVn);
ARR22=-inv(inv(AWh)*AWp+inv(AVh)*AVp)*(inv(AWh)*AWn+inv(AVh)*AVn);

