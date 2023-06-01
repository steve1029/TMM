
alpha       = sqrt(ep.*mu - ep(1)*mu(1)*(sin(inc_ang))^2);
xsi         = exp(i*k0*alpha.*d);
if mode == 0% TE
    eta     = alpha./mu;
else        % TM
    eta     = alpha./ep;
end

%%
A           = zeros(1,N);
B           = zeros(1,N);
A(1)        = 1;


b=zeros(1,N);
inv_a=zeros(1,N);

%% at output
F=1;
G=eta(N);

%% from last to first
for q=N-1:-1:2
    temp=inv([1 1;eta(q) -eta(q)]) * [F;G];
    b(q)=temp(2);
    inv_a(q)=1/temp(1);
    xsiBA1xsi=xsi(q)*b(q)*inv_a(q)*xsi(q);
    F=   1   * (1+xsiBA1xsi);
    G=eta(q) * (1-xsiBA1xsi);    
end

%% at input
temp=inv([F -1;G eta(1)]) * [1;eta(1)];
A(2)        = temp(1);
B(1)        = temp(2);

%% from first to last
for q=2:N-1
    A(q+1)  = inv_a(q)*xsi(q)*A(q);
    B(q)    = b(q)*A(q+1);
end


%%

trans=A(N);
refle=B(1);
