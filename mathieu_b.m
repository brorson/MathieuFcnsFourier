function b = mathieu_b(m, q)
  % This returns the Mathieu eigenvalue ("characteristic value") 
  % for order m and frequency parameter q.

  N = m+25; % max(10,2*m);
    
  tol = 1e-14;
  if (abs(mod(m,2)) < tol)
    % Even order m
    MM = make_matrix_oe(N,q);
    [V,D] = eig(MM);  
    [D,idx] = sort(diag(D));

    col = round(m/2);
    %fprintf('mathieu_b order %d -- extracting col = %d\n', m, col)
    b = D(col);
  else
    % Odd order m
    MM = make_matrix_oo(N,q);
    [V,D] = eig(MM);  
    [D,idx] = sort(diag(D));

    col = round((m-1)/2+1);
    %fprintf('mathieu_b order %d -- extracting col = %d\n', m, col)
    b = D(col);

  end
  
