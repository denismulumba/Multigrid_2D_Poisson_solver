function v_f = prolong(nx, ny, v)
    % interpolate 'v' to the fine grid
    v_f = zeros(2*nx+2, 2*ny+2);
 for i = 2 : nx
     for j = 2: ny
         v_f(2*i-1,2*j-1) = 0.5625*v(i,j)+0.1875*(v(i-1,j)+v(i,j-1))+0.0625*v(i-1,j-1);
         v_f(2*i  ,2*j-1) = 0.5625*v(i,j)+0.1875*(v(i+1,j)+v(i,j-1))+0.0625*v(i+1,j-1);
         v_f(2*i-1,2*j  ) = 0.5625*v(i,j)+0.1875*(v(i-1,j)+v(i,j+1))+0.0625*v(i-1,j+1);
         v_f(2*i  ,2*j  ) = 0.5625*v(i,j)+0.1875*(v(i+1,j)+v(i,j+1))+0.0625*v(i+1,j+1);
     end
 end

end
