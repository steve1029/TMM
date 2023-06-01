% case 1. homogeneous space - finite PC waveguide - homogeneous space 

SNlay=10;  %  the period number of finite PC waveguide

Ta=zeros(2*L,2*L,SNlay); % left to rignt
Ra=zeros(2*L,2*L,SNlay); % left to right

Tb=zeros(2*L,2*L,SNlay); % right to left
Rb=zeros(2*L,2*L,SNlay); % right to left

Ca=zeros(4*L,2*L,SNlay); % left to right
Cb=zeros(4*L,2*L,SNlay); % left to right


tCa=zeros(4*L,2*L,SNlay); % left to right
tCb=zeros(4*L,2*L,SNlay); % left to right

for laynt=1:SNlay
    
    Ra(:,:,laynt)=AR11;
    Ta(:,:,laynt)=AT12;
    Rb(:,:,laynt)=AR22;
    Tb(:,:,laynt)=AT21;
    
    Ca(:,:,laynt)=ACa;
    Cb(:,:,laynt)=ACb;
    
end;


   I=eye(2*L,2*L);
   T_temp1a=Ta(:,:,1);
   R_temp1a=Ra(:,:,1);
   T_temp1b=Tb(:,:,1);
   R_temp1b=Rb(:,:,1);
   
   %%% Important
   tCa=Ca;
   tCb=Cb;
   %%%

for laynt=2:SNlay
   
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
   
   SAT12=T_temp1a;
   SAR11=R_temp1a;
   SAT21=T_temp1b;
   SAR22=R_temp1b;
   SACa=Ca;
   SACb=Cb;
   
%% Field visualization

LFMM_Amode_field_visual_xz_case1_Lfree_Rfree_leftright;
% LFMM_Amode_field_visual_xz_case1_Lfree_Rfree_rightleft;





   
   
   
   
   