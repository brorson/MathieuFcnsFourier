function A = mathieu_coeffs_ee(N,q,m)
  % Returns Fourier coeffs for the mth order ce_2n Mathieu fcn.
  % Allowed value of m = 0, 2, 4, 6, ... 

  tol = 1e-14;
  if (abs(mod(m,2)) > tol)
    error('Must invoke ee function with even order m')
  end
  
  M = make_matrix_ee(N+4,q);
  
  %fprintf('mathieu_coeffs_ee, condeig(M) = \n')
  %disp(condeig(M))
  
  % This invokes arpack
  %[V,D] = eigs(M,N,'largestabs');

  % This invokes lapack
  [V,D] = eig(M);  

  % Sort ascending so low-order A coeffs are on LHS of
  % V matrix.
  [D,idx] = sort(diag(D));
  V = V(:,idx);

  %fprintf('D(1) = %f, D(2) = %f\n', D(1), D(2))  
  
  % Undo change to matrix made in make_matrix_ee
  V(1,:) = V(1,:)/sqrt(2);


  col = round(m/2+1);
  s = 1; %sign(V(1,col));
  %fprintf('mathieu_coeffs_ee, col = %d\n', col)
  A = s*V(1:N,col);
  
end

