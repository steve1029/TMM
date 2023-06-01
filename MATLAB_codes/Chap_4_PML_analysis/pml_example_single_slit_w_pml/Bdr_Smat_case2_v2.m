% boundary matching S-matrix
% 2. homogeneous space         - grating - inhomogeneous waveguide

%% input + body
inv2a1b=inv_func(RRa,Lfree_Rf);
inv1b2a=inv_func(Lfree_Rf,RRa);
for k=1:Nlay;[Ca(:,:,k) Cb(:,:,k)]  =Redheffer_Cr(Lfree_Rb,Lfree_Tf,Lfree_Rf,Lfree_Tb,RRa,TTa,RRb,TTb,inv2a1b,inv1b2a,Ca(:,:,k),Cb(:,:,k));end;
[RRa TTa RRb TTb]                   =Redheffer_RT(Lfree_Rb,Lfree_Tf,Lfree_Rf,Lfree_Tb,RRa,TTa,RRb,TTb,inv2a1b,inv1b2a);
 
%% body + output
inv2a1b=inv_func(Rwg_Rb2,RRb);
inv1b2a=inv_func(RRb,Rwg_Rb2);
for k=1:Nlay;[Ca(:,:,k) Cb(:,:,k)]  =Redheffer_Cl(RRa,TTa,RRb,TTb,Rwg_Rb2,Rwg_Tf2,Rwg_Rf2,Rwg_Tb2,inv2a1b,inv1b2a,Ca(:,:,k),Cb(:,:,k));end;
[RRa TTa RRb TTb]                   =Redheffer_RT(RRa,TTa,RRb,TTb,Rwg_Rb2,Rwg_Tf2,Rwg_Rf2,Rwg_Tb2,inv2a1b,inv1b2a);

