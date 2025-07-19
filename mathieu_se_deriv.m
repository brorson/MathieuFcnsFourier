function sed = mathieu_se_deriv(m, q, v)
  % This computes the derivative of the Mathieu se 
  % fcn order m for frequency parameter q and angle v.
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
    %fprintf('Even Mathieu se deriv, m = %d\n', m)    
    B = mathieu_coeffs_oe(N,q,m);
    s = ones(size(B));
    % sum on 2, 4, 6, 8 ...
    % Sum from smallest to largest coeff
    % k = (1:N)';
    sedm = zeros(size(v));
    sedp = zeros(size(v));
    for k=N:-1:1
      if (s(k)<0)
	sedm = sedm + 2*k*B(k)*cos(2*k*v);
      else
	sedp = sedp + 2*k*B(k)*cos(2*k*v);	
      end
      sed = sedp+sedm;
    end
    % Hack -- make sure fcn is positive at v = 0.
    %if (q<0)
    %  if (mod(m,4) < tol)
    %    s(2:2:end) = -1;
    %  else
    %    s(1:2:end) = -1;
    %  end
    %end
    %s = sign(sum(s.*B));
    %sed = s*sed;
 
  else
    % Odd
    %fprintf('Odd Mathieu se deriv, m = %d\n', m)
    B = mathieu_coeffs_oo(N,q,m);
    s = ones(size(B));
    % sum on 1, 3, 5, 7, ...
    % Sum from smallest to largest coeff
    % k = (0:N-1)';
    sedp = zeros(size(v));
    sedm = zeros(size(v));    
    for k=(N-1):-1:0
      if (s(k+1)<0)
	sedm = sedm + (2*k+1)*B(k+1)*cos((2*k+1)*v);
      else
	sedp = sedp + (2*k+1)*B(k+1)*cos((2*k+1)*v);	
      end
    end
    sed = sedm+sedp;
    % Hack -- make sure fcn has the right sign.
    if (q<0)
      if (mod(m-1,4) < tol)
        s(2:2:end) = -1;
      else
        s(1:2:end) = -1;
      end
    end
    s = sign(sum(s.*B));
    sed = s*sed;
  end

  % By defintion, all se slopes are positive for v = 0.  Modify fcns
  % to obey this definition.
  %if (B(1) < 0)
  %  sed = -sed;
  %end
  
end
