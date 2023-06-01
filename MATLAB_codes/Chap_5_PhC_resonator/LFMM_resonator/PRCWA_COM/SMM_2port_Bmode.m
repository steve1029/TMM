% SMM_2port_Bmode

Iden=zeros(2*L,2*L);
for k=1:2*L
Iden(k,k)=1;   
end;

KII=zeros(2*L,2*L); % KII matrix
for k=1:NBx
   for l=1:NBy
od_ind1=(k-1)*NBy+l;      
K11(od_ind1,od_ind1)=-(Bkz_vc(k)*Bky_vc(l))/Bkx_vc(k,l);
K12(od_ind1,od_ind1)=-(Bkx_vc(k,l)^2+Bkz_vc(k)^2)/Bkx_vc(k,l);
K21(od_ind1,od_ind1)=(Bky_vc(l)^2+Bkx_vc(k,l)^2)/Bkx_vc(k,l);
K22(od_ind1,od_ind1)=( Bky_vc(l)*Bkz_vc(k) )/Bkx_vc(k,l);
   end;
end;

KII(1:L, 1:L)=K11;
KII(1:L, L+1:2*L)=K12;
KII(L+1:2*L, 1:L)=K21;
KII(L+1:2*L,L+1:2*L )=K22;

BYh=Iden;
BZh=KII*(1/(w0*mu0));

BYpm=zeros(2*L,pcnt); BZpm=zeros(2*L,pcnt); 
BYnm=zeros(2*L,mcnt); BZnm=zeros(2*L,mcnt);

BYpp=zeros(2*L,pcnt); BZpp=zeros(2*L,pcnt);
BYnp=zeros(2*L,mcnt); BZnp=zeros(2*L,mcnt);


for pbm_ind=1:pcnt % BYpm, BZpm, BYpp, BZpp

    for k=1:NBx
        for l=1:NBy
   
    [ey ez]=byp(xm-xm,k,l,pbm_ind,BpEFy,BpEFz,bp_evalue);
    BYpm(NBy*(k-1)+l,pbm_ind)=ey;
    BYpm(L+NBy*(k-1)+l,pbm_ind)=ez;
    
    [hy hz]=bzp(xm-xm,k,l,pbm_ind,BpHFy,BpHFz,bp_evalue);
    BZpm(NBy*(k-1)+l,pbm_ind)=hy;
    BZpm(L+NBy*(k-1)+l,pbm_ind)=hz;
    
    [ey ez]=byp(xp-xm,k,l,pbm_ind,BpEFy,BpEFz,bp_evalue);
    BYpp(NBy*(k-1)+l,pbm_ind)=ey;
    BYpp(L+NBy*(k-1)+l,pbm_ind)=ez;
    
    [hy hz]=bzp(xp-xm,k,l,pbm_ind,BpHFy,BpHFz,bp_evalue);
    BZpp(NBy*(k-1)+l,pbm_ind)=hy;
    BZpp(L+NBy*(k-1)+l,pbm_ind)=hz;
    
        end;
    end;

end;

for pbm_ind=1:mcnt % AWnm, AVnm, AWnp, AVnp
    
    for k=1:NBx
        for l=1:NBy
   
    [ey ez]=byn(xm-xp,k,l,pbm_ind,BnEFy,BnEFz,bm_evalue);
    BYnm(NBy*(k-1)+l,pbm_ind)=ey;
    BYnm(L+NBy*(k-1)+l,pbm_ind)=ez;
    
    [hy hz]=bzn(xm-xp,k,l,pbm_ind,BnHFy,BnHFz,bm_evalue);
    BZnm(NBy*(k-1)+l,pbm_ind)=hy;
    BZnm(L+NBy*(k-1)+l,pbm_ind)=hz;
    
    [ey ez]=byn(xp-xp,k,l,pbm_ind,BnEFy,BnEFz,bm_evalue);
    BYnp(NBy*(k-1)+l,pbm_ind)=ey;
    BYnp(L+NBy*(k-1)+l,pbm_ind)=ez;
    
    [hy hz]=bzn(xp-xp,k,l,pbm_ind,BnHFy,BnHFz,bm_evalue);
    BZnp(NBy*(k-1)+l,pbm_ind)=hy;
    BZnp(L+NBy*(k-1)+l,pbm_ind)=hz;
        end;
    end;

end;

% Layer S-matrix 
U=eye(2*NBx*NBy);
S11=inv(BYh)*BYpm+inv(BZh)*BZpm;
S12=inv(BYh)*BYnm+inv(BZh)*BZnm;
S21=inv(BYh)*BYpp-inv(BZh)*BZpp;
S22=inv(BYh)*BYnp-inv(BZh)*BZnp;
S=[S11 S12;S21 S22];

%left-to-right layer S-matrix
BCa=inv(S)*[2*U;zeros(2*NBx*NBy)];  %% very good!!!
BCap=BCa(1:pcnt,:);
BCam=BCa(pcnt+1:pcnt+mcnt,:);
BR33=inv(BYh)*(BYpm*BCap+BYnm*BCam-BYh*U);
BT34=inv(BYh)*(BYpp*BCap+BYnp*BCam);

% right-to-left layer S-matrix
BCb=inv(S)*[zeros(2*NBx*NBy);2*U];
BCbp=BCb(1:pcnt,:);
BCbm=BCb(pcnt+1:pcnt+mcnt,:);
BR44=inv(BYh)*(BYpp*BCbp+BYnp*BCbm-BYh*U);
BT43=inv(BYh)*(BYpm*BCbp+BYnm*BCbm);

% Boundary S-matrix

BYp=zeros(2*L,pcnt); BZp=zeros(2*L,pcnt); 
BYn=zeros(2*L,mcnt); BZn=zeros(2*L,mcnt);

for pbm_ind=1:pcnt  % BYpm, BZpm, BYpp, BZpp

    for k=1:NBx
        for l=1:NBy
   
    [ey ez]=byp(xc,k,l,pbm_ind,BpEFy,BpEFz,bp_evalue);
    BYp(NBy*(k-1)+l,pbm_ind)=ey;
    BYp(L+NBy*(k-1)+l,pbm_ind)=ez;
    
    [hy hz]=bzp(xc,k,l,pbm_ind,BpHFy,BpHFz,bp_evalue);
    BZp(NBy*(k-1)+l,pbm_ind)=hy;
    BZp(L+NBy*(k-1)+l,pbm_ind)=hz;

    
        end;
    end;

end;

for pbm_ind=1:mcnt % BYnm, BZnm, BYnp, BZnp
    
    for k=1:NBx
        for l=1:NBy
   
    [ey ez]=byn(xc,k,l,pbm_ind,BnEFy,BnEFz,bm_evalue);
    BYn(NBy*(k-1)+l,pbm_ind)=ey;
    BYn(L+NBy*(k-1)+l,pbm_ind)=ez;
    
    [hy hz]=bzn(xc,k,l,pbm_ind,BnHFy,BnHFz,bm_evalue);
    BZn(NBy*(k-1)+l,pbm_ind)=hy;
    BZn(L+NBy*(k-1)+l,pbm_ind)=hz;

    
        end;
    end;

end;

% left side semi-infinite waveguide
BLR33=-inv(inv(BYh)*BYn-inv(BZh)*BZn)*(inv(BYh)*BYp-inv(BZh)*BZp);
BLT34=inv(inv(BYn)*BYh-inv(BZn)*BZh)*(inv(BYn)*BYp-inv(BZn)*BZp);
BLT43=2*inv(inv(BYh)*BYn-inv(BZh)*BZn);
BLR44=-inv(inv(BYn)*BYh-inv(BZn)*BZh)*(inv(BYn)*BYh+inv(BZn)*BZh);

% right side semi-infinite waveguide
BRR33=-inv(inv(BYp)*BYh+inv(BZp)*BZh)*(inv(BYp)*BYh-inv(BZp)*BZh);
BRT34=2*inv(inv(BYh)*BYp+inv(BZh)*BZp);
BRT43=inv(inv(BYp)*BYh+inv(BZp)*BZh)*(inv(BYp)*BYn-inv(BZp)*BZn);
BRR44=-inv(inv(BYh)*BYp+inv(BZh)*BZp)*(inv(BYh)*BYn+inv(BZh)*BZn);


