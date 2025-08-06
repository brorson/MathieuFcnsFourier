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

    sep = zeros(size(v));
    sem = zeros(size(v));
    sedm = zeros(size(v));
    sedp = zeros(size(v));
    
    % sum on 2, 4, 6, 8 ... sum pos and neg terms separately for stability
    for k=N:-1:1
      t = B(k)*sin(2*k*v);
      td = 2*k*B(k)*cos(2*k*v);
      if (t<0)
        sem = sem + t;
      else
        sep = sep + t;
      end
      if (td<0)
        sedm = sedm + td;
      else
        sedp = sedp + td;
      end
    end
    % I should do a sort before doing the sum
    se = sep+sem;
    sed = sedp+sedm;
    % Hack -- make sure deriv is positive at v = 0.
    s = ones(size(B));
    %if (q<0)
    %  if (mod(m-1,4) < tol)
    %    s(2:2:end) = -1;
    %  else
    %    s(1:2:end) = -1;
    %  end
    %end
    ss = sign(sum(s.*B));
    se = ss.*se;
    sed = ss.*sed;
    
  else
    % Odd
    % fprintf('Odd Mathieu se, m = %d\n', m)
    B = mathieu_coeffs_oo(N,(q),m);

    sep = zeros(size(v));
    sem = zeros(size(v));
    sedp = zeros(size(v));
    sedm = zeros(size(v));    
    % sum on 1, 3, 5, 7, ... sum pos and neg terms separately for stability
    for k=(N-1):-1:0
      t = B(k+1)*sin((2*k+1)*v);
      td = (2*k+1)*B(k+1)*cos((2*k+1)*v);
      if (t<0)
        sem = sem + t;
      else
        sep = sep + t;
      end
    
      if (td<0)
        sedm = sedm + td;
      else
        sedp = sedp + td;
      end
    end
    se = sep+sem;
    sed = sedp+sedm;
    
    % Hack -- make sure deriv is positive at v = 0 for q<0.
    s = ones(size(B));
    if (q<0)
      if (mod(m-1,4) < tol)
        s(2:2:end) = -1;
      else
        s(1:2:end) = -1;
      end
    end
    ss = sign(sum(s.*B));
    se = ss.*se;
    sed = ss.*sed;
  end

end
