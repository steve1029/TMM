%% space determiniation
[m,n]                   = meshgrid([1:1:num_hx]-NBx,[1:1:num_hy]-NBy);
m                       = m';
n                       = n';

eps_L                   = zeros(num_hx,num_hy,NN);
aps_L                   = zeros(num_hx,num_hy,NN);
mu_L                    = zeros(num_hx,num_hy,NN);
bu_L                    = zeros(num_hx,num_hy,NN);

full_pattern            = rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);

%% pattern for the first half
for nn = 1 : NN/2

    xn                  = a/NN * nn;

    hole_pa             = zeros(num_hx,num_hy);

    if xn < rr
        p_size          = sqrt(rr^2 - xn^2);
        for pp = 1 : PP
            offset      = sqrt(3)*a*(pp-1) + sqrt(3)*a*first_offset/2;
            hole_pa     = hole_pa + ...
                rect_2D_mesh(m,n,1,Tx,Ty,+offset-p_size,+offset+p_size,-Ty/2,Ty/2) + ...
                rect_2D_mesh(m,n,1,Tx,Ty,-offset-p_size,-offset+p_size,-Ty/2,Ty/2);
        end
        offset          = sqrt(3)*a*((PP+1)-1) + sqrt(3)*a*first_offset/2; % should be equal to Tx/2
        hole_pa         = hole_pa + ...
            rect_2D_mesh(m,n,1,Tx,Ty,+offset-p_size,+offset,-Ty/2,Ty/2) + ...
            rect_2D_mesh(m,n,1,Tx,Ty,-offset,-offset+p_size,-Ty/2,Ty/2);
    end
    
    if abs(xn - a/2) < rr
        q_size          = sqrt(rr^2 - (xn - a/2)^2);
        for qq = 1 : QQ
            offset      = sqrt(3)*a*(qq-1) + sqrt(3)*a*first_offset/2 + sqrt(3)*a/2;
            hole_pa    = hole_pa + ...
                rect_2D_mesh(m,n,1,Tx,Ty,+offset-q_size,+offset+q_size,-Ty/2,Ty/2) + ...
                rect_2D_mesh(m,n,1,Tx,Ty,-offset-q_size,-offset+q_size,-Ty/2,Ty/2);
        end
%         hole_pa         = hole_pa + ... % to test the complete bandgap
%             rect_2D_mesh(m,n,1,Tx,Ty,-q_size,+q_size,-Ty/2,Ty/2);
    end
    
    inverse_pa          = full_pattern - hole_pa;
    
    eps_L(:,:,nn)       =   epra * hole_pa +   eprd * inverse_pa;
    aps_L(:,:,nn)       = 1/epra * hole_pa + 1/eprd * inverse_pa;
    mu_L (:,:,nn)       =   mur0 * full_pattern;
    bu_L (:,:,nn)       = 1/mur0 * full_pattern;    
end

%% copy to the other half side
for nn = NN : -1 : NN/2+1
    eps_L(:,:,nn)       = eps_L(:,:,NN+1-nn);
    aps_L(:,:,nn)       = aps_L(:,:,NN+1-nn);
    mu_L(:,:,nn)        = mu_L (:,:,NN+1-nn);
    bu_L(:,:,nn)        = bu_L (:,:,NN+1-nn);

end


