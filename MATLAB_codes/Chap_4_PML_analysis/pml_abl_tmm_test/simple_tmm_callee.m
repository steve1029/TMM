
alpha       = sqrt(ep.*mu - ep(1)*mu(1)*(sin(inc_ang))^2);
xsi         = exp(i*k0*alpha.*d);
if mode == 0% TE
    eta     = alpha./mu;
else        % TM
    eta     = alpha./ep;
end

M_L         = zeros(2,2,N);
M_R         = zeros(2,2,N);

for q=1:N
    M_L(:,:,q)   = [xsi(q)       1;  xsi(q)*eta(q)         -eta(q)];
    M_R(:,:,q)   = [     1  xsi(q);         eta(q)  -xsi(q)*eta(q)];
end

temp        = eye(2);
for q=1:N-1
    temp    = temp * inv(M_L(:,:,q)) * M_R(:,:,q+1); % µÚ¿¡ q+1 !!
end

trans       = 1        /temp(1,1);
refle       = temp(2,1)/temp(1,1);


