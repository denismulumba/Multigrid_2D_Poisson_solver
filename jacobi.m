

function [u_new, u] = jacobi(u, f, h, omega, Npre, Npost)
% Performs Npre pre-smoothing and Npost post-smoothing Jacobi iterations on u

% Initialize solution guess
    u_new = u;
    nlevels = 8;
    NX = 1*2^(nlevels-1); % Nx and Ny are given as function of grid levels
    NY = 1*2^(nlevels-1);
    % Pre-smoothing
    for i = 1:Npre
        u_old = u_new;
        for j = 2:size(u,1)-1
            for k = 2:size(u,2)-1
                u_new(j,k) = (1-omega)*u_old(j,k) + omega*(0.25*(u_old(j-1,k) + u_old(j+1,k) + u_old(j,k-1) + u_old(j,k+1)) - h^2*f(j,k));
            end
        end
    end
    
    % Post-smoothing
    for i = 1:Npost
        u_old = u_new;
        for j = 2:size(u,1)-1
            for k = 2:size(u,2)-1
                u_new(j,k) = (1-omega)*u_old(j,k) + omega*(0.25*(u_old(j-1,k) + u_old(j+1,k) + u_old(j,k-1) + u_old(j,k+1)) - h^2*f(j,k));
            end
        end
    end
    res = zeros(NX+2, NY+2);
    for i=2:NX+1
        for j=2:NY+1
            res(i,j) = f(i,j) - ((Ax*(u(i+1,j)+u(i-1,j)) + Ay*(u(i,j+1)+u(i,j-1)) - 2.0*(Ax+Ay)*u(i,j)));
        end
    end

end