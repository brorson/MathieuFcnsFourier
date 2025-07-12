function a = mathieu_a(m, q)
  % This returns the Mathieu eigenvalue ("characteristic value") 
  % for order m and   frequency parameter q.

  N = max(10,2*m);
    
  tol = 1e-14;
  if (abs(mod(m,2)) < tol)
    % Even order m
    MM = make_matrix_ee(N,q);
    [V,D] = eig(MM);  
    [D,idx] = sort(diag(D));

    col = round(m/2+1);
    %fprintf('mathieu_a order %d -- extracting col = %d\n', m, col)
    a = D(col);
  else
    % Odd order m
    MM = make_matrix_eo(N,q);
    [V,D] = eig(MM);  
    [D,idx] = sort(diag(D));

    col = round((m+1)/2);
    %fprintf('mathieu_a order %d -- extracting col = %d\n', m, col)
    a = D(col);

  end
  
