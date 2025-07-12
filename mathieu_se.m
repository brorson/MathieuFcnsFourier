function se = mathieu_se(m, q, v)
  % This computes the Mathieu se fcn order m for frequency
  % parameter q and angle v.
  % Input m, q must be scalars, v may be a vector.
    
  % I find the peak Fourier coeff tracks m.  Therefore
  % Adjust the matrix size based on order m.
  N = m+10;
  
  % Use different coeffs depending upon whether m is even or
  % odd.
  tol = 1e-14;
  if (abs(mod(m,2) < tol))
    % Even
    fprintf('Even Mathieu se, m = %d\n', m)    
    B = mathieu_coeffs_oe(N,q,m);
    % sum on 2, 4, 6, 8 ...
    % Sum from smallest to largest coeff
    % k = (1:N)';
    se = zeros(size(v));
    for k=N:-1:1
      se = se + B(k)*sin(2*k*v);
    end
  else
    % Odd
    fprintf('Odd Mathieu se, m = %d\n', m)
    B = mathieu_coeffs_oo(N,q,m);
    % sum on 1, 3, 5, 7, ...
    % Sum from smallest to largest coeff
    % k = (0:N-1)';
    se = zeros(size(v));
    for k=(N-1):-1:0
      se = se + B(k+1)*sin((2*k+1)*v);
    end
  end

  % By defintion, all slopes are positive for v = 0.  Modify fcns
  % to obey this definition.
  if (B(1) < 0)
    se = -se;
  end
  
end
