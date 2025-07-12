function ced = mathieu_ce_deriv(m, q, v)
  % This computes the derivative of the Mathieu fcn ce 
  % order m for frequency parameter q and angle v.
  % The derivative is taken in the v (angular) direction.
  % Input m, q must be scalars, v may be a vector.
    
  % I find the peak Fourier coeff tracks m.  Therefore
  % Adjust the matrix size based on order m.
  N = m+10;

  % Use different coeffs depending upon whether m is even or
  % odd.
  tol = 1e-14;
  if (abs(mod(m,2) < tol))
    % Even
    fprintf('Even Mathieu ce deriv, m = %d\n', m)    

    A = mathieu_coeffs_ee(N,q,m);
    % Sum on 0, 2, 4, 6, ...
    % Sum from smallest to largest coeff.
    % k = (0:N-1)';
    ced = zeros(size(v));
    for k=(N-1):-1:0
      ced = ced - 2*k*A(k+1)*sin(2*k*v);
    end
  else
    % Odd
    fprintf('Odd Mathieu ce deriv, m = %d\n', m)
    A = mathieu_coeffs_eo(N,q,m);
    % Sum on 1, 3, 5, 7, ...
    % Sum from smallest to largest coeff.
    % k = (0:N-1)';
    ce = zeros(size(v));
    for k=(N-1):-1:0
      ced = ced - (2*k+1)*A(k+1)*sin((2*k+1)*v);
    end
  end

  % By defintion, all fcns are positive for v = 0.  Modify fcn
  % to obey this definition.
  if (A(1) < 0)
    ced = -ced;
  end
    
end
