% boundary matching S-matrix
% 1. homogeneous space - grating - homogeneous space       

% Final boundary matching
%% input + body
    I=eye(2*L,2*L);
    T_temp1a=Lfree_Tf;       %outTb
    R_temp1a=Lfree_Rb;       %outRf
    T_temp1b=Lfree_Tb;       %outTf
    R_temp1b=Lfree_Rf;       %outRb
 
 %%% Important
    tCa=Ca;
    tCb=Cb;
 %%%
    %    
    T_temp2a=TTa;
    R_temp2a=RRa;
    T_temp2b=TTb;
    R_temp2b=RRb;
 
  RRa=(R_temp1a+T_temp1b*inv(I-R_temp2a*R_temp1b)*R_temp2a*T_temp1a);
 	TTa=T_temp2a*inv(I-R_temp1b*R_temp2a)*T_temp1a;
% 
 	RRb=(R_temp2b+T_temp2a*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b);
    TTb=T_temp1b*inv(I-R_temp2a*R_temp1b)*T_temp2b;
   
    for k=1:Nlay
    tCa(:,:,k)=Ca(:,:,k)*inv(I-R_temp1b*R_temp2a)*T_temp1a;
    tCb(:,:,k)=Cb(:,:,k)+Ca(:,:,k)*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b;
  	end;
    Ca=tCa;
    Cb=tCb;
 
 %% body + output
    T_temp1a=TTa;
    R_temp1a=RRa;
    T_temp1b=TTb;
    R_temp1b=RRb;
%    
    T_temp2a=Rfree_Tf;
    R_temp2a=Rfree_Rb;
    T_temp2b=Rfree_Tb;
    R_temp2b=Rfree_Rf;
    
    RRa=(R_temp1a+T_temp1b*inv(I-R_temp2a*R_temp1b)*R_temp2a*T_temp1a);
 	TTa=T_temp2a*inv(I-R_temp1b*R_temp2a)*T_temp1a;
% 
 	RRb=(R_temp2b+T_temp2a*inv(I-R_temp1b*R_temp2a)*R_temp1b*T_temp2b);
    TTb=T_temp1b*inv(I-R_temp2a*R_temp1b)*T_temp2b;
   
   
    for k=1:Nlay
%    	
    tCa(:,:,k)=Ca(:,:,k)+Cb(:,:,k)*inv(I-R_temp2a*R_temp1b)*R_temp2a*T_temp1a;
    tCb(:,:,k)=Cb(:,:,k)*inv(I-R_temp2a*R_temp1b)*T_temp2b;
% 
    end; % for k
    Ca=tCa;
    Cb=tCb;

 
%% part 5 Field visualization & data save

switch direct_
    
    case 1      % left-to-right
        
        Field_visualization_3D_xz_case1_Lfree_Rfree_leftright;

    case 2      % right-to-left
        
        Field_visualization_3D_xz_case1_Lfree_Rfree_rightleft;

end;















