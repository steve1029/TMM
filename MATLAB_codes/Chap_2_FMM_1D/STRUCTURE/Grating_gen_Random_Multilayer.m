        
   eps_L=zeros(1,Nlay);
   aps_L=zeros(1,Nlay);
   mu_L =zeros(1,Nlay);
   bu_L =zeros(1,Nlay);
  
   eps_L(1,:)=(1+rand(1,Nlay)+0.0000001*i).^2; % random permittivity
   aps_L=1./eps_L;
   mu_L(1,:)=1;  % random permeability
   bu_L=1./mu_L;
   
   eps_L(1,1)   =(ni+0.000001*i)^2;
   eps_L(1,Nlay)=(nf+0.000001*i)^2;
  
