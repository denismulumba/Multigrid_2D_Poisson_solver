

function v_c = restrict(nx, ny, v)
% This function restricts a fine grid solution to a coarse grid using the bilinear interpolation scheme
% nx and ny are the dimensions of the coarse grid
% v is the fine grid solution to be restricted

% Initialize the coarse grid solution with zeros
v_c=zeros([nx+2,ny+2]);

% Interpolate the fine grid solution onto the coarse grid using bilinear interpolation
for i = 2: nx+1
    for j = 2: ny+1
        v_c(i,j)=0.25*(v(2*i-1,2*j-1)+v(2*i,2*j-1)+v(2*i-1,2*j)+v(2*i,2*j));
    end
end
end
