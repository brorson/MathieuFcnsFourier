function B = mathieu_coeffs_oe(N,q,m)
  % Returns Fourier coeffs for the mth order se_2n Mathieu fcn.
  % Allowed value of m = 2, 4, 6, ... 

  tol = 1e-14;
  if (abs(mod(m,2)) > tol)
    error('Must invoke oe function with even order m')
  end
  
  M = make_matrix_oe(2*N,q);
  
  %[V,D] = eigs(M,N,'largestabs');
  [V,D] = eig(M);  
  
  [D,idx] = sort(diag(D));
  V = V(:,idx);
  %disp(D(idx))

  col = round(m/2);
  fprintf('mathieu_coeffs_oe, col = %d\n', col)
  B = V(1:N,col);
  
end

