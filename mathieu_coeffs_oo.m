function B = mathieu_coeffs_oo(N,q,m)
  % Returns Fourier coeffs for the mth order se_2n+1 Mathieu fcn.
  % Allowed value of m = 1, 3, 5, ... 

  tol = 1e-14;
  if (abs(mod(m+1,2)) > tol)
    error('Must invoke oo function with odd order m')
  end
  
  M = make_matrix_oo(2*N,q);
  
  %[V,D] = eigs(M,N,'largestabs');
  [V,D] = eig(M);  
  
  [D,idx] = sort(diag(D));
  V = V(:,idx);
  %disp(D(idx))

  col = round((m-1)/2+1);
  s = 1;% sign(V(1,col));
  %fprintf('mathieu_coeffs_oo, col = %d\n', col)
  B = s*V(1:N,col);
  
end

