function A = mathieu_coeffs_ee(N,q,m)
  % Returns Fourier coeffs for the mth order ce_2n Mathieu fcn.
  % Allowed value of m = 0, 2, 4, 6, ... 

  tol = 1e-14;
  if (abs(mod(m,2)) > tol)
    error('Must invoke ee function with even order m')
  end
  
  M = make_matrix_ee(2*N,q);
  
  %[V,D] = eigs(M,N,'largestabs');
  [V,D] = eig(M);  
  
  [D,idx] = sort(diag(D));
  V = V(:,idx);
  %disp(D(idx))
  
  % Undo change to matrix
  V(1,:) = V(1,:)/sqrt(2);

  col = round(m/2+1);
  fprintf('mathieu_coeffs_ee, col = %d\n', col)
  A = V(1:N,col);
  
end

