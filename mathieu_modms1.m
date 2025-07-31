function [Ms,Msd] = mathieu_modms1(m, q, u)
  % This computes the modified odd Mathieu fcn of the 
  % first kind of order m for frequency
  % parameter q and radius u.  This is the fcn frequently
  % called Ms in the literature.
  % Inputs m, q must be scalars, u may be a vector.
    
  % Here I use the expression given in Zhang and Jin.

  if (q<0)
    error('Modified Mathieu fcns for negatative q not implemented yet!\n')
  end

  if (m<1)
    error('Invalid order m requested for Ms1!\n')
  end

  % Force v to be col vector so returns are col vectors.
  if (size(u, 2)>1)
    u = u';
  end

    
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
    %fprintf('Even Mathieu Ms(1), m = %d\n', m)
    B = mathieu_coeffs_oe(N,q,m);
    Ms = 0;
    Msd = 0;
    for k=(N-1):-1:0
      Jks = besselj(k,s);
      Jkp2t = besselj(k+2,t);
      Jkp2s = besselj(k+2,s);
      Jkt = besselj(k,t);

      Jdks = besseljd(k,s);
      Jdkp2t = besseljd(k+2,t);
      Jdkp2s = besseljd(k+2,s);
      Jdkt = besseljd(k,t);

      Ms = Ms + ((-1).^k).*B(k+1).*(Jks.*Jkp2t - Jkp2s.*Jkt);
      Msd = Msd + ((-1)^k)*B(k+1)*...
	    (exp(u).*(Jks.*Jdkp2t - Jkp2s.*Jdkt) - ...
	     exp(-u).*(Jdks.*Jkp2t - Jdkp2s.*Jkt));
    end
    % Do normalization.
    sgn = (m-2)/2;
    Ms = (((-1)^sgn)/B(1))*Ms;    
    Msd = sqrt(q)*(((-1)^sgn)/B(1))*Msd; 
  else
    % Odd
    %fprintf('Odd Mathieu Ms(1), m = %d\n', m)
    B = mathieu_coeffs_oo(N,q,m);
    Ms = 0;
    Msd = 0;    
    for k=(N-1):-1:0    
      Jks = besselj(k,s);
      Jkt = besselj(k,t);
      Jkp1s = besselj(k+1,s);
      Jkp1t = besselj(k+1,t);

      Jdks = besseljd(k,s);
      Jdkt = besseljd(k,t);
      Jdkp1s = besseljd(k+1,s);
      Jdkp1t = besseljd(k+1,t);

      Ms = Ms + ((-1).^k).*B(k+1).*(Jks.*Jkp1t - Jkp1s.*Jkt);
      Msd = Msd + ((-1)^k)*B(k+1)*...
	    (exp(u).*(Jks.*Jdkp1t - Jkp1s.*Jdkt) - ...
	     exp(-u).*(Jdks.*Jkp1t - Jdkp1s.*Jkt));
    end
    % Do normalization.
    sgn = (m-1)/2;
    Ms = (((-1)^sgn)/B(1))*Ms;    
    Msd = sqrt(q)*(((-1)^sgn)/B(1))*Msd; 
  end

  % At this point, fcns should be properly normalized.
  
end
