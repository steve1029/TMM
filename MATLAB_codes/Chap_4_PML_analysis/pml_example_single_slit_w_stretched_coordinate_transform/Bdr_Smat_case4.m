% boundary matching S-matrix
% 4. inhomogeneous waveguide  - grating - inhomogeneous waveguide

% Final boundary matching
 % input + body
    I=eye(2*L,2*L);
    T_temp1a=Lwg_Tf1;       %outTb
    R_temp1a=Lwg_Rb1;       %outRf
    T_temp1b=Lwg_Tb1;       %outTf
    R_temp1b=Lwg_Rf1;       %outRb
 
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
    T_temp2a=Rwg_Tf2;
    R_temp2a=Rwg_Rb2;
    T_temp2b=Rwg_Tb2;
    R_temp2b=Rwg_Rf2;
    
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

%Field_visualization_3D_xz_case4_Lwg_Rwg_leftright;
%Field_visualization_3D_xz_case4_Lwg_Rwg_rightleft;

