function Ms = mathieu_modms1(m, q, u)
  % This computes the modified odd Mathieu fcn of the 
  % first kind of order m for frequency
  % parameter q and radius u.  This is the fcn frequently
  % called Ms in the literature.
  % Inputs m, q must be scalars, u may be a vector.
    
  % Here I use the expression given in Briscombe, Corless,
  % et al., https://arxiv.org/abs/2008.01812, eq. 5.7.

  % I find the peak Fourier coeff tracks m.  Therefore
  % I adjust the matrix size based on order m.
  N = m+10;

  % Utility vars.
  sqq = sqrt(q);
  s = sqq*exp(u);
  t = sqq*exp(-u);
  
  % Use different coeffs depending upon whether m is even or
  % odd.
  tol = 1e-14;
  if (abs(mod(m,2) < tol))
    % Even
    fprintf('Even Mathieu Ms, m = %d\n', m)
    B = mathieu_coeffs_oe(N,q,m);
    Ms = 0;
    for k=N:-1:1
      Jkp1s = besselj(k+1,s);
      Jkm1t = besselj(k-1,t);
      Jkm1s = besselj(k-1,s);
      Jkp1t = besselj(k+1,t);
      Ms = Ms + ((-1).^k).*B(k).*(Jkp1s.*Jkm1t - Jkm1s.*Jkp1t);
    end
    % Do normalization.
    %Ms = (((-1)^m)/B(1))*sqrt(pi/2)*Ms;
    Ms = (((-1)^m)/B(1))*Ms;    

  else
    % Odd
    fprintf('Odd Mathieu Ms, m = %d\n', m)
    B = mathieu_coeffs_eo(N,q,m);
    Ms = 0;
    for k=(N-1):-1:0    
      Jks = besselj(k,s);
      Jkt = besselj(k,t);
      Jkp1s = besselj(k+1,s);
      Jkp1t = besselj(k+1,t);
      Ms = Ms + ((-1).^k).*B(k+1).*(Jkp1s.*Jkt - Jks.*Jkp1t);
    end
    % Do normalization.
    %Ms = (((-1)^m)/B(1))*sqrt(pi/2)*Ms;
    Ms = (((-1)^m)/B(1))*Ms;    
  end

  % By defintion, all derivs are positive for u = 0.  Modify fcn
  % to obey this definition.
  %[~, zidx] = min(abs(u));
  %if ((Ms(zidx+1)-Ms(zidx))<0)
  %  Ms = -Ms;
  %end

  % At this point, fcns should be properly normalized.
  
end
