function Mc = mathieu_modmc2(m, q, u)
  % This computes the modified Mathieu fcn of the 
  % second kind of order m for frequency
  % parameter q and radius u.  This is the fcn frequently
  % called Mc in the literature.
  % Inputs m, q must be scalars, u may be a vector.
    
  % Here I use the expressions given in Zhang and Jin.

  % I find the peak Fourier coeff tracks m.  Therefore
  % I adjust the matrix size based on order m.
  N = m+10;

  % Utility vars.
  sqq = sqrt(q);
  s = sqq*exp(-u);
  t = sqq*exp(u);

  % Use different coeffs depending upon whether m is even or
  % odd.
  tol = 1e-14;
  if (abs(mod(m,2) < tol))
    % Even
    fprintf('Even Mathieu Mc, m = %d\n', m)
    A = mathieu_coeffs_ee(N,q,m);
    Mc = 0;
    for k=(N-1):-1:0
      Jks = besselj(k,s);
      Ykt = bessely(k,t);
      Mc = Mc + ((-1)^k)*A(k+1)*Jks.*Ykt ;
    end
    % Do normalization.
    sgn = m/2;
    Mc = (((-1)^sgn)/A(1))*Mc;    
  else
    % Odd
    fprintf('Odd Mathieu Mc, m = %d\n', m)
    A = mathieu_coeffs_eo(N,q,m);
    Mc = 0;
    for k=(N-1):-1:0
      Jks = besselj(k,s);
      Ykt = bessely(k,t);
      Jkp1s = besselj(k+1,s);
      Ykp1t = bessely(k+1,t);
      Mc = Mc + ((-1)^k)*A(k+1)*(Jks.*Ykp1t + Jkp1s.*Ykt);
    end
    % Do normalization.
    sgn = (m-1)/2;
    Mc = (((-1)^sgn)/A(1))*Mc;    
  end

  % By defintion, all fcns are positive for u = 0.  Modify fcn
  % to obey this definition.
  %[~, zidx] = min(abs(u));
  %if (Mc(zidx)<0)
  %  Mc = -Mc;
  %end
  %[~,k] = max(abs(A));
  %if (A(k)*((-1)^(k-1)) < 0)
  %  %Mc = -Mc;
  %end

  % At this point, fcns should be properly normalized.
  
end
