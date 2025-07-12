function A = make_matrix_eo(N, q)
  % This creates the recurrence relation matrix for
  % the even-odd Mathieu fcns (ce_2n+1).
  % Inputs:
  % N = matrix size (related to max order desired).
  % q = shape parameter.
  % Output:
  % A = sparse recurrence matrix.

  if (N < 3)
    error('Must ask for matrix of size 3 or greater.')
  end
  
  % Form matrix
  %A = sparse(N,N);  
  A = zeros(N,N);    
  A(1,1) = 1+q;
  A(1,2) = q;
  A(2,1) = q;
  A(2,2) = 3^2;
  A(2,3) = q;
  
  for j = 3:N-1
    A(j,j-1) = q;
    A(j,j)   = (2*j-1)^2;
    A(j,j+1) = q;
  end
  
  A(N,N-1) = q;
  A(N,N)   = (2*N-1)^2;

end

