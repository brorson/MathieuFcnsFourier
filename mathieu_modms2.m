function Ms = mathieu_modms2(m, q, u)
  % This computes the modified odd Mathieu fcn of the 
  % second kind of order m for frequency
  % parameter q and radius u.  This is the fcn frequently
  % called Ms in the literature.
  % Inputs m, q must be scalars, u may be a vector.
    
  % Here I use the expression given in Zhang and Jin.

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
    fprintf('Even Mathieu Ms, m = %d\n', m)
    B = mathieu_coeffs_oe(N,q,m);
    Ms = 0;
    for k=(N-1):-1:0
      Jks = besselj(k,s);
      Ykt = bessely(k,t);
      Jkp2s = besselj(k+2,s);
      Ykp2t = bessely(k+2,t);
      Ms = Ms + ((-1).^k).*B(k+1).*(Jks.*Ykp2t - Jkp2s.*Ykt);
    end
    % Do normalization.
    sgn = (m-2)/2;
    Ms = (((-1)^sgn)/B(1))*Ms;    

  else
    % Odd
    fprintf('Odd Mathieu Ms, m = %d\n', m)
    B = mathieu_coeffs_oo(N,q,m);
    Ms = 0;
    for k=(N-1):-1:0    
      Jks = besselj(k,s);
      Ykt = bessely(k,t);
      Jkp1s = besselj(k+1,s);
      Ykp1t = bessely(k+1,t);
      Ms = Ms + ((-1).^k).*B(k+1).*(Jks.*Ykp1t - Jkp1s.*Ykt);
    end
    % Do normalization.
    sgn = (m-1)/2;
    Ms = (((-1)^sgn)/B(1))*Ms;    
  end

  % By defintion, all derivs are positive for u = 0.  Modify fcn
  % to obey this definition.
  %[~, zidx] = min(abs(u));
  %if ((Ms(zidx+1)-Ms(zidx))<0)
  %  Ms = -Ms;
  %end

  % At this point, fcns should be properly normalized.
  
end
