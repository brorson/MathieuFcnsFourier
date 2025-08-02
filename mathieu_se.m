function [se, sed] = mathieu_se(m, q, v)
  % This computes the Mathieu se fcn order m for frequency
  % parameter q and angle v.
  % Input m, q must be scalars, v may be a vector.
  % Returns:  se = Mathieu se fcn, sed = deriv of se fcn.
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
    %fprintf('Even Mathieu se, m = %d\n', m)    

    B = mathieu_coeffs_oe(N,(q),m);
    s = ones(size(B));
    % sum on 2, 4, 6, 8 ... sum pos and neg terms separately for stability
    sep = zeros(size(v));
    sem = zeros(size(v));
    sedm = zeros(size(v));
    sedp = zeros(size(v));
    for k=N:-1:1
      if (s(k)<0)
        sem = sem + s(k)*B(k)*sin(2*k*v);
	sedm = sedm + 2*k*B(k)*cos(2*k*v);
      else
        sep = sep + s(k)*B(k)*sin(2*k*v);
	sedp = sedp + 2*k*B(k)*cos(2*k*v);
      end
    end
    se = sep-sem;
    sed = sedp-sedm;
    % Hack -- make sure deriv is positive at v = 0.
    s = sign(sum(s.*B));
    se = s*se;
    sed = s*sed;
    
  else
    % Odd
    % fprintf('Odd Mathieu se, m = %d\n', m)
    B = mathieu_coeffs_oo(N,(q),m);
    s = ones(size(B));
    % sum on 1, 3, 5, 7, ... sum pos and neg terms separately for stability
    sep = zeros(size(v));
    sem = zeros(size(v));
    sedp = zeros(size(v));
    sedm = zeros(size(v));    
    for k=(N-1):-1:0
      if (s(k+1)<0)
        sem = sem + s(k+1)*B(k+1)*sin((2*k+1)*v);
	sedm = sedm + (2*k+1)*B(k+1)*cos((2*k+1)*v);
      else
        sep = sep + s(k+1)*B(k+1)*sin((2*k+1)*v);
	sedp = sedp + (2*k+1)*B(k+1)*cos((2*k+1)*v);
      end
    end
    se = sep-sem;
    sed = sedp-sedm;
    % Hack -- make sure deriv is positive at v = 0 for q<0.
    if (q<0)
      if (mod(m-1,4) < tol)
        s(2:2:end) = -1;
      else
        s(1:2:end) = -1;
      end
    end
    s = sign(sum(s.*B));
    se = s*se;
    sed = s*sed;
  end

end
