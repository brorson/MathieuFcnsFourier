function A = mathieu_coeffs_eo(N,q,m)
  % Returns Fourier coeffs for the mth order ce_2n+1 Mathieu fcn.
  % Allowed value of m = 1, 3, 5, 7, ... 

  tol = 1e-14;
  if (abs(mod(m+1,2)) > tol)
    error('Must invoke eo function with odd order m')
  end
  
  M = make_matrix_eo(N,q);
  
  %fprintf('mathieu_coeffs_eo, condeig(M) = \n')
  %disp(condeig(M))

  %[V,D] = eigs(M,N,'largestabs');
  [V,D] = eig(M);  
  
  [D,idx] = sort(diag(D));
  V = V(:,idx);
  %disp(D(idx))
  
  col = round( (m+1)/2 );
  %fprintf('mathieu_coeffs_eo, col = %d\n', col)
  A = V(1:N,col);
  
end
