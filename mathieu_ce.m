function [ce, ced] = mathieu_ce(m, q, v)
  % This computes the Mathieu fcn ce order m for frequency
  % parameter q and angle v.
  % Input m, q must be scalars, v may be a vector.
  % Returns:  ce = Mathieu ce fcn, ced = deriv of ce fcn.
  % The returns are col vectors.
  
  % Force v to be col vector so returns are col vectors.
  if (size(v, 2)>1)
    v = v';
  end
  
    
  % I find the peak Fourier coeff tracks m.  Therefore
  % Adjust the matrix size based on order m.
  N = m+25;

  % Use different coeffs depending upon whether m is even or
  % odd.
  tol = 1e-14;
  if (abs(mod(m,2) < tol))
    % Even
    %fprintf('Even Mathieu ce, m = %d\n', m)    

    A = mathieu_coeffs_ee(N,(q),m);
    s = ones(size(A));
    % Sum on 0, 2, 4, 6, ...
    % Sum from smallest to largest coeff.
    % k = (0:N-1)';
    cem = zeros(size(v));  % Actual fcn -- minus terms
    cep = zeros(size(v));  % Acutal fcn -- pos terms
    cedm = zeros(size(v)); % Deriv -- minus terms
    cedp = zeros(size(v)); % Deriv -- plus terms
    for k=(N-1):-1:0
      if (s(k+1)<0)
        cem = cem + s(k+1)*A(k+1)*cos(2*k*v);
	cedm = cedm - 2*k*A(k+1)*sin(2*k*v);
      else
        cep = cep + s(k+1)*A(k+1)*cos(2*k*v);
	cedp = cedp - 2*k*A(k+1)*sin(2*k*v);  
      end
    end
    % I should do a sort before doing the sum
    ce = cep+cem;
    ced = cedp+cedm;
    % Hack -- this makes sure the fcn has the right overall sign.
    if (q<0)
      if (mod(m,4) < tol)
        s(2:2:end) = -1;
      else
        s(1:2:end) = -1;
      end
    end
    s = sign(sum(s.*A));
    ce = s*ce;
    ced = s*ced;

  else
    % Odd
    %fprintf('Odd Mathieu ce, m = %d\n', m)

    A = mathieu_coeffs_eo(N,(q),m);
    s = ones(size(A));
    % Sum on 1, 3, 5, 7, ...
    % Sum from smallest to largest coeff.
    % k = (0:N-1)';
    cem = zeros(size(v));
    cep = zeros(size(v));
    cedm = zeros(size(v));
    cedp = zeros(size(v));
    for k=(N-1):-1:0
      if (s(k+1)<0)
        cem = cem + s(k+1)*A(k+1)*cos((2*k+1)*v);
	cedm = cedm - (2*k+1)*A(k+1)*sin((2*k+1)*v);
      else
        cep = cep + s(k+1)*A(k+1)*cos((2*k+1)*v);
	cedp = cedp - (2*k+1)*A(k+1)*sin((2*k+1)*v);
      end
    end
    ce = cep+cem;
    ced = cedp+cedm;
    % Hack -- make sure fcn has correct overall sign.
    if (q<0)
      if (mod(m-1,4) < tol)
        s(2:2:end) = -1;
      else
        s(1:2:end) = -1;
      end
    end
    s = sign(sum(s.*A));
    ce = s*ce;
    ced = s*ced;
  end

end
