function [Ms,Msd] = mathieu_modms2(m, q, u)
  % This computes the modified odd Mathieu fcn of the 
  % second kind of order m for frequency
  % parameter q and radius u.  This is the fcn frequently
  % called Ms in the literature.
  % Inputs m, q must be scalars, u may be a vector.
    
  % Here I use the expression given in Zhang and Jin.

  if (q<0)
    error('Modified Mathieu fcns for negatative q not implemented yet!\n')
  end

  if (m<1)
    error('Invalid order m requested for Ms2!\n')
  end

  % Force v to be col vector so returns are col vectors.
  if (size(u, 2)>1)
    u = u';
  end

 
  % Set offset used in Bessel fcn depending upon order m.
  % This is per the book "Accurately Calculating Mathieu Functions",
  % XXXXX & YYYY
  c = 2;
  
    
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
    %fprintf('Even Mathieu Ms(2), m = %d\n', m)
    B = mathieu_coeffs_oe(N,q,m);
    Ms = 0;
    Msd = 0;
    for k=(N-1):-1:0
      Jks = besselj(k,s);
      Ykt = bessely(k,t);
      Jkp2s = besselj(k+2,s);
      Ykp2t = bessely(k+2,t);

      Jdks = besseljd(k,s);
      Ydkp2t = besselyd(k+2,t);
      Jdkp2s = besseljd(k+2,s);
      Ydkt = besselyd(k,t);
      
      Ms = Ms + ((-1).^k).*B(k+1).*(Jks.*Ykp2t - Jkp2s.*Ykt);
      Msd = Msd + ((-1)^k)*B(k+1)*...
            (exp(u).*(Jks.*Ydkp2t - Jkp2s.*Ydkt) - ...
             exp(-u).*(Jdks.*Ykp2t - Jdkp2s.*Ykt));
    end
    % Do normalization.
    sgn = (m-2)/2;
    Ms = (((-1)^sgn)/B(1))*Ms;    
    Msd = sqrt(q)*(((-1)^sgn)/B(1))*Msd; 
  else
    % Odd
    %fprintf('Odd Mathieu Ms(2), m = %d\n', m)
    B = mathieu_coeffs_oo(N,q,m);
    Ms = 0;
    Msd = 0;
    for k=(N-1):-1:0    
      Jks = besselj(k,s);
      Ykt = bessely(k,t);
      Jkp1s = besselj(k+1,s);
      Ykp1t = bessely(k+1,t);

      Jdks = besseljd(k,s);
      Ydkt = besselyd(k,t);
      Jdkp1s = besseljd(k+1,s);
      Ydkp1t = besselyd(k+1,t);
      
      Ms = Ms + ((-1).^k).*B(k+1).*(Jks.*Ykp1t - Jkp1s.*Ykt);
      Msd = Msd + ((-1)^k)*B(k+1)*...
            (exp(u).*(Jks.*Ydkp1t - Jkp1s.*Ydkt) - ...
             exp(-u).*(Jdks.*Ykp1t - Jdkp1s.*Ykt));
    end
    % Do normalization.
    sgn = (m-1)/2;
    Ms = (((-1)^sgn)/B(1))*Ms;
    Msd = sqrt(q)*(((-1)^sgn)/B(1))*Msd;     
  end

  % By defintion, all derivs are positive for u = 0.  Modify fcn
  % to obey this definition.
  %[~, zidx] = min(abs(u));
  %if ((Ms(zidx+1)-Ms(zidx))<0)
  %  Ms = -Ms;
  %end

  % At this point, fcns should be properly normalized.
  
end
