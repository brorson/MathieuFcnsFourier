function Ce = mathieu_modce1(m, q, u)
  % This computes the modified Mathieu fcn of the 
  % first kind of order m for frequency
  % parameter q and radius u.  This is the fcn frequently
  % called Ce in the literature.
  % Inputs m, q must be scalars, u may be a vector.
    
  % Here I use the expressions given in Briscombe, Corless,
  % et al., https://arxiv.org/abs/2008.01812

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
    fprintf('Even Mathieu Ce, m = %d\n', m)
    A = mathieu_coeffs_ee(N,q,m);
    Ce = 0;
    for k=(N-1):-1:0
      Jks = besselj(k,s);
      Jkt = besselj(k,t);
      Ce = Ce + ((-1)^k)*A(k+1)*Jks.*Jkt ;
    end
    % Do normalization.
    %Ce = (((-1)^m)/A(1))*sqrt(pi/2)*Ce;
    Ce = (((-1)^m)/A(1))*Ce;    
  else
    % Odd
    fprintf('Odd Mathieu Ce, m = %d\n', m)
    A = mathieu_coeffs_eo(N,q,m);
    Ce = 0;
    for k=(N-1):-1:0
      Jks = besselj(k,s);
      Jkt = besselj(k,t);
      Jkp1s = besselj(k+1,s);
      Jkp1t = besselj(k+1,t);
      Ce = Ce + ((-1)^k)*A(k+1)*(Jkp1s.*Jkt + Jks.*Jkp1t);
    end
    % Do normalization.
    % Ce = (((-1)^m)/A(1))*sqrt(pi/2)*Ce;
    Ce = (((-1)^m)/A(1))*Ce;    
  end

  % By defintion, all fcns are positive for u = 0.  Modify fcn
  % to obey this definition.
  %[~, zidx] = min(abs(u));
  %if (Ce(zidx)<0)
  %  Ce = -Ce;
  %end
  %[~,k] = max(abs(A));
  %if (A(k)*(-1^(k-1)) < 0)
  %  %Ce = -Ce;
  %end

  % At this point, fcns should be properly normalized.
  
end
