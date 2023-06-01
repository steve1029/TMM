%% space determiniation
[m,n]                   = meshgrid([1:1:num_hx]-NBx,[1:1:num_hy]-NBy);
m                       = m';
n                       = n';

eps_L                   = zeros(num_hx,num_hy,Nlay);
aps_L                   = zeros(num_hx,num_hy,Nlay);
mu_L                    = zeros(num_hx,num_hy,Nlay);
bu_L                    = zeros(num_hx,num_hy,Nlay);

full_pattern            = rect_2D_mesh(m,n,1,Tx,Ty,-Tx/2,Tx/2,-Ty/2,Ty/2);

%% pattern for the first half
offset_C                = 8 * nm;
offset_B                = 16 * nm;
offset_A                = 24 * nm;

for nn = 1 : Nlay/2

    xn                  = cell_num*a/Nlay * nn;

    hole_pa             = zeros(num_hx,num_hy);

    if abs(xn - 0*a/2) < rr                             % region (1)
        p_size          = sqrt(rr^2 - (xn - 0*a/2)^2);
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
    
    if abs(xn - 1*a/2) < rr                             % region (2)
        q_size          = sqrt(rr^2 - (xn - 1*a/2)^2);
        for qq = 1 : QQ
            offset      = sqrt(3)*a*(qq-1) + sqrt(3)*a*first_offset/2 + sqrt(3)*a/2;
            hole_pa    = hole_pa + ...
                rect_2D_mesh(m,n,1,Tx,Ty,+offset-q_size,+offset+q_size,-Ty/2,Ty/2) + ...
                rect_2D_mesh(m,n,1,Tx,Ty,-offset-q_size,-offset+q_size,-Ty/2,Ty/2);
        end
    end
    
    if abs(xn - 2*a/2) < rr                             % region (3)
        p_size          = sqrt(rr^2 - (xn - 2*a/2)^2);
        offset          = sqrt(3)*a*( 1-1) + sqrt(3)*a*first_offset/2 + offset_C;
        hole_pa         = hole_pa + ...
            rect_2D_mesh(m,n,1,Tx,Ty,+offset-p_size,+offset+p_size,-Ty/2,Ty/2) + ...
            rect_2D_mesh(m,n,1,Tx,Ty,-offset-p_size,-offset+p_size,-Ty/2,Ty/2);
        for pp = 2 : PP
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

    if abs(xn - 3*a/2) < rr                             % region (4)
        q_size          = sqrt(rr^2 - (xn - 3*a/2)^2);
        offset          = sqrt(3)*a*( 1-1) + sqrt(3)*a*first_offset/2 + sqrt(3)*a/2 + offset_C;
        hole_pa         = hole_pa + ...
            rect_2D_mesh(m,n,1,Tx,Ty,+offset-q_size,+offset+q_size,-Ty/2,Ty/2) + ...
            rect_2D_mesh(m,n,1,Tx,Ty,-offset-q_size,-offset+q_size,-Ty/2,Ty/2);
        for qq = 2 : QQ
            offset      = sqrt(3)*a*(qq-1) + sqrt(3)*a*first_offset/2 + sqrt(3)*a/2;
            hole_pa    = hole_pa + ...
                rect_2D_mesh(m,n,1,Tx,Ty,+offset-q_size,+offset+q_size,-Ty/2,Ty/2) + ...
                rect_2D_mesh(m,n,1,Tx,Ty,-offset-q_size,-offset+q_size,-Ty/2,Ty/2);
        end
    end

    if abs(xn - 4*a/2) < rr                             % region (5)
        p_size          = sqrt(rr^2 - (xn - 4*a/2)^2);
        offset          = sqrt(3)*a*( 1-1) + sqrt(3)*a*first_offset/2 + offset_B;
        hole_pa         = hole_pa + ...
            rect_2D_mesh(m,n,1,Tx,Ty,+offset-p_size,+offset+p_size,-Ty/2,Ty/2) + ...
            rect_2D_mesh(m,n,1,Tx,Ty,-offset-p_size,-offset+p_size,-Ty/2,Ty/2);
        offset          = sqrt(3)*a*( 2-1) + sqrt(3)*a*first_offset/2 + offset_C;
        hole_pa         = hole_pa + ...
            rect_2D_mesh(m,n,1,Tx,Ty,+offset-p_size,+offset+p_size,-Ty/2,Ty/2) + ...
            rect_2D_mesh(m,n,1,Tx,Ty,-offset-p_size,-offset+p_size,-Ty/2,Ty/2);
        for pp = 3 : PP
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

    if abs(xn - 5*a/2) < rr                             % region (6)
        q_size          = sqrt(rr^2 - (xn - 5*a/2)^2);
        offset          = sqrt(3)*a*( 1-1) + sqrt(3)*a*first_offset/2 + sqrt(3)*a/2 + offset_B;
        hole_pa         = hole_pa + ...
            rect_2D_mesh(m,n,1,Tx,Ty,+offset-q_size,+offset+q_size,-Ty/2,Ty/2) + ...
            rect_2D_mesh(m,n,1,Tx,Ty,-offset-q_size,-offset+q_size,-Ty/2,Ty/2);
        for qq = 2 : QQ
            offset      = sqrt(3)*a*(qq-1) + sqrt(3)*a*first_offset/2 + sqrt(3)*a/2;
            hole_pa    = hole_pa + ...
                rect_2D_mesh(m,n,1,Tx,Ty,+offset-q_size,+offset+q_size,-Ty/2,Ty/2) + ...
                rect_2D_mesh(m,n,1,Tx,Ty,-offset-q_size,-offset+q_size,-Ty/2,Ty/2);
        end
    end

    if abs(xn - 6*a/2) < rr                             % region (7)
        p_size          = sqrt(rr^2 - (xn - 6*a/2)^2);
        offset          = sqrt(3)*a*( 1-1) + sqrt(3)*a*first_offset/2 + offset_A;
        hole_pa         = hole_pa + ...
            rect_2D_mesh(m,n,1,Tx,Ty,+offset-p_size,+offset+p_size,-Ty/2,Ty/2) + ...
            rect_2D_mesh(m,n,1,Tx,Ty,-offset-p_size,-offset+p_size,-Ty/2,Ty/2);
        offset          = sqrt(3)*a*( 2-1) + sqrt(3)*a*first_offset/2 + offset_C;
        hole_pa         = hole_pa + ...
            rect_2D_mesh(m,n,1,Tx,Ty,+offset-p_size,+offset+p_size,-Ty/2,Ty/2) + ...
            rect_2D_mesh(m,n,1,Tx,Ty,-offset-p_size,-offset+p_size,-Ty/2,Ty/2);
        for pp = 3 : PP
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

    if abs(xn - 7*a/2) < rr                             % region (8)
        q_size          = sqrt(rr^2 - (xn - 7*a/2)^2);
        offset          = sqrt(3)*a*( 1-1) + sqrt(3)*a*first_offset/2 + sqrt(3)*a/2 + offset_B;
        hole_pa         = hole_pa + ...
            rect_2D_mesh(m,n,1,Tx,Ty,+offset-q_size,+offset+q_size,-Ty/2,Ty/2) + ...
            rect_2D_mesh(m,n,1,Tx,Ty,-offset-q_size,-offset+q_size,-Ty/2,Ty/2);
        for qq = 2 : QQ
            offset      = sqrt(3)*a*(qq-1) + sqrt(3)*a*first_offset/2 + sqrt(3)*a/2;
            hole_pa    = hole_pa + ...
                rect_2D_mesh(m,n,1,Tx,Ty,+offset-q_size,+offset+q_size,-Ty/2,Ty/2) + ...
                rect_2D_mesh(m,n,1,Tx,Ty,-offset-q_size,-offset+q_size,-Ty/2,Ty/2);
        end
    end

    inverse_pa          = full_pattern - hole_pa;
    
    eps_L(:,:,nn)       =   epra * hole_pa +   eprd * inverse_pa;
    aps_L(:,:,nn)       = 1/epra * hole_pa + 1/eprd * inverse_pa;
    mu_L (:,:,nn)       =   mur0 * full_pattern;
    bu_L (:,:,nn)       = 1/mur0 * full_pattern;    
end

%% copy to the other half side
for nn = Nlay : -1 : Nlay/2+1
    eps_L(:,:,nn)       = eps_L(:,:,Nlay+1-nn);
    aps_L(:,:,nn)       = aps_L(:,:,Nlay+1-nn);
    mu_L(:,:,nn)        = mu_L (:,:,Nlay+1-nn);
    bu_L(:,:,nn)        = bu_L (:,:,Nlay+1-nn);

end


