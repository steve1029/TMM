%% RCWA field visualization x-y plane

xx=[-Tx/2:Tx*0.01:Tx/2];
yy=[-Ty/2:Ty*0.01:Ty/2];

[x,y]=meshgrid(xx,yy);

Ex_xy=zeros(length(xx),length(yy));
Ey_xy=zeros(length(xx),length(yy));
Ez_xy=zeros(length(xx),length(yy));
Hx_xy=zeros(length(xx),length(yy));
Hy_xy=zeros(length(xx),length(yy));
Hz_xy=zeros(length(xx),length(yy));

 zv=2*um;
 
    for k=1:NBx
      for l=1:NBy
         
        
         Ex_xy=Ex_xy+ET3x_cof(k,l)*exp(j*( kz_vc(k,l)*(zv-ac_thick(Nlay)) ))*exp(j*( kx_vc(k)*x+ky_vc(l)*y ));
         Ey_xy=Ey_xy+ET3y_cof(k,l)*exp(j*( kz_vc(k,l)*(zv-ac_thick(Nlay)) ))*exp(j*( kx_vc(k)*x+ky_vc(l)*y ));
         Ez_xy=Ez_xy+ET3z_cof(k,l)*exp(j*( kz_vc(k,l)*(zv-ac_thick(Nlay)) ))*exp(j*( kx_vc(k)*x+ky_vc(l)*y ));
         
         Hx_xy=Hx_xy+HT3x_cof(k,l)*exp(j*( kz_vc(k,l)*(zv-ac_thick(Nlay)) ))*exp(j*( kx_vc(k)*x+ky_vc(l)*y ));
         Hy_xy=Hy_xy+HT3y_cof(k,l)*exp(j*( kz_vc(k,l)*(zv-ac_thick(Nlay)) ))*exp(j*( kx_vc(k)*x+ky_vc(l)*y ));
         Hz_xy=Hz_xy+HT3z_cof(k,l)*exp(j*( kz_vc(k,l)*(zv-ac_thick(Nlay)) ))*exp(j*( kx_vc(k)*x+ky_vc(l)*y ));
         
         
      end;
	end;

% evanescent field distribution
figure(1); mesh(abs((Ex_xy)));
figure(2); mesh(abs((Ey_xy)));
figure(3); mesh(abs((Ez_xy)));

% propagating field distribution
figure(5); mesh(abs((Hx_xy)));
figure(6); mesh(abs((Hy_xy)));
figure(7); mesh(abs((Hz_xy)));

 
%  
% if zv <= 0 % region I
%    
%    Eincx=zeros(length(xx),length(yy));
%    Eincy=zeros(length(xx),length(yy));
%    Eincz=zeros(length(xx),length(yy));
%    
%    phase_inc = exp(i*(kix*x + kiy*y + kiz*zv));%
%    
%    Eincx = Ux*phase_inc ;
%    Eincy = Uy*phase_inc ;
%    Eincz = Uz*phase_inc ;
% 
%    for k=1:NBx
%       for l=1:NBy
%          total_Ex_xy=total_Ex_xy+Rx_cof(k,l)*exp(j*( kx_ref(k)*x+ky_ref(l)*y-kz_ref(k,l)*zv ));
%          total_Ey_xy=total_Ey_xy+Ry_cof(k,l)*exp(j*( kx_ref(k)*x+ky_ref(l)*y-kz_ref(k,l)*zv ));
% 		 total_Ez_xy=total_Ez_xy+Rz_cof(k,l)*exp(j*( kx_ref(k)*x+ky_ref(l)*y-kz_ref(k,l)*zv ));
% 		end;
% 	end;
%    
%    total_Ex_xy=total_Ex_xy+Eincx;
%    total_Ey_xy=total_Ey_xy+Eincy;
%    total_Ez_xy=total_Ez_xy+Eincz;
% 
%    
% elseif  zv < ac_thick(Nlay) % grating region
% 
% 
% 		for ni=1:Nlay
% 				if (zv>=ac_thick(ni)-lay_thick(ni))&(zv<ac_thick(ni))
%    			vi=ni;
%    			end % if
% 		end;
%       
%             vi
% 
% 			X=zeros(2*L,2*L);	
% 			Xinv=zeros(2*L,2*L);
%             Xtemp=zeros(4*L,4*L);
%             for k=1:2*L
% 			X(k,k)=exp( Teigvalue(1,k,vi)*(zv-(ac_thick(vi)-lay_thick(vi))) );   
%             Xinv(k,k)=exp( -Teigvalue(1,k,vi)*(zv-(ac_thick(vi))) );  
%          	end;
%          	Xtemp(1:2*L,1:2*L)=X;
%             Xtemp(2*L+1:4*L,2*L+1:4*L)=Xinv;
%             
%             H=zeros(4*L,4*L);
% 			H(1:2*L,1:2*L)=WW(:,:,vi);
%             H(1:2*L,2*L+1:4*L)=WW(:,:,vi);
% 			H(2*L+1:4*L,1:2*L)=j*k0*VV(:,:,vi);
% 			H(2*L+1:4*L,2*L+1:4*L)=-j*k0*VV(:,:,vi);
% 
% 
%             SUtemp1=H*Xtemp;
%             SUtemp2=SUtemp1*Lay_cof(:,vi);
%             tzy=SUtemp2(2*L+1:3*L,1);
%             tzx=SUtemp2(3*L+1:4*L,1);
%             SUz=-1*inv(Ez)*(Kx*tzy-Ky*tzx)/k0;
% 
%             
%             for k=1:NBx
%                for l=1:NBy
%             total_Ey_xy=total_Ey_xy+SUtemp2((k-1)*N_SNx+l,1)*exp(j*( kx_ref(k)*x+ky_ref(l)*y )); %Sx
%             total_Ex_xy=total_Ex_xy+SUtemp2(L+(k-1)*N_SNx+l,1)*exp(j*( kx_ref(k)*x+ky_ref(l)*y )); %Sy
%             total_Ez_xy=total_Ez_xy+SUz((k-1)*N_SNy+l)*exp(j*( kx_ref(k)*x+ky_ref(l)*y ));
%          		end; % for l
%       		end; % for k
% 
%       
% else % region II
%    
%    for k=1:NBx
%       for l=1:NBy
%          total_Ex_xy=total_Ex_xy+Tx_cof(k,l)*exp(j*( kx_tra(k)*x+ky_tra(l)*y+kz_ref(k,l)*(zv-ac_thick(Nlay)) ));
%          total_Ey_xy=total_Ey_xy+Ty_cof(k,l)*exp(j*( kx_tra(k)*x+ky_tra(l)*y+kz_ref(k,l)*(zv-ac_thick(Nlay)) ));
% 		 total_Ez_xy=total_Ez_xy+Tz_cof(k,l)*exp(j*( kx_tra(k)*x+ky_tra(l)*y+kz_ref(k,l)*(zv-ac_thick(Nlay)) ));
% 		end;
% 	end;
% 
% 
% 	        
% end; % if
% 
% figure(1); imagesc(abs(real(total_Ey_xy)));
% figure(2); imagesc(abs(real(total_Ex_xy)));
% figure(3); imagesc(abs(real(total_Ez_xy)));

