

function [u, res] = V_cycle(nx, ny, num_levels, u, f, level)
% V cycle
if nargin < 6
level = 1;
end

if (level == num_levels) % bottom solve
[u, res] = Jacrelax(level, nx, ny, u, f, 1, true);
return
end

% Step 1: Relax Au=f on this grid
[u, res] = Jacrelax(level, nx, ny, u, f, 1, true);

% Step 2: Restrict residual to coarse grid
res_c = restrict(floor(nx/2), floor(ny/2), res);

% Step 3: Solve A e_c=res_c on the coarse grid.
e_c = zeros(size(res_c));
[e_c, res_c] = V_cycle(nx/2, floor(ny/2), num_levels, e_c, res_c, level+1); 

% Step 4: Interpolate(prolong) e_c to fine grid and add to u
u = u + prolong(floor(nx/2), floor(ny/2), e_c);

% Step 5: Relax Au=f on this grid 
[u, res] = Jacrelax(level, nx, ny, u, f, 1, false);

end




