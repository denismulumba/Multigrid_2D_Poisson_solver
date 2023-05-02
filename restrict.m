function v_c = restrict(nx, ny, v)
  % restrict 'v' to the coarser grid
%   v_c = zeros(nx+2, ny+2);
%   v_c(2:nx+1, 2:ny+1) = 0.25*(v(1:2:end-1,1:2:end-1) + v(2:2:end,1:2:end-1) + ...
%  v(1:2:end-1,2:2:end)   + v(2:2:end,2:2:end));

  v_c=zeros([nx+2,ny+2]);

%vectorized form of
  for i = 2: nx+1
      for j = 2: ny+1
      v_c(i,j)=0.25*(v(2*i-1,2*j-1)+v(2*i,2*j-1)+v(2*i-1,2*j)+v(2*i,2*j));

      end
  end
end

