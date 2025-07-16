function A = make_matrix_oe(N, q)
  % This creates the recurrence relation matrix for
  % the odd-even Mathieu fcns (se_2n).
  % Inputs:
  % N = matrix size (related to max order desired).
  % q = shape parameter.
  % Output:
  % A = recurrence matrix.
    
  if (N < 3)
    error('Must ask for matrix of size 3 or greater.')
  end
    
  % Form matrix
  %A = sparse(N,N);  
  A = zeros(N,N);    
  A(1,1) = 2^2;
  A(1,2) = q;
  A(2,1) = q;
  A(2,2) = 4^2;
  A(2,3) = q;
  
  for j = 3:N-1
    A(j,j-1) = q;
    A(j,j)   = (2*j)^2;
    A(j,j+1) = q;
  end
  
  A(N,N-1) = q;
  A(N,N)   = (2*N)^2;

end

