function [Mc,Mcd] = mathieu_modmc1(m, q, u)
  % This computes the modified Mathieu fcn of the 
  % first kind of order m for frequency
  % parameter q and radius u.  This is the fcn frequently
  % called Mc in the literature.  It also computes its
  % u derivative.
  % Inputs m, q must be scalars, u may be a vector.
    
  % Here I use the expressions given in Zhang and Jin

  if (q<0)
    error('Modified Mathieu fcns for negatative q not implemented yet!\n')
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
    % Even -- m = 0, 2, 4, 6, ...
    %fprintf('Even Mathieu Mc(1), m = %d\n', m)
    A = mathieu_coeffs_ee(N,q,m);
    Mc = 0;
    Mcd = 0;
    for k=(N-1):-1:0
      Jks = besselj(k,s);
      Jkt = besselj(k,t);
      Jdks = besseljd(k,s);
      Jdkt = besseljd(k,t);
      Mc = Mc + ((-1)^k)*A(k+1)*Jks.*Jkt ;
      Mcd = Mcd + ((-1)^k)*A(k+1)*(exp(u).*Jks.*Jdkt - exp(-u).*Jdks.*Jkt) ;
    end
    % Do normalization.
    sgn = m/2;
    Mc = (((-1)^sgn)/A(1))*Mc;    
    Mcd = sqrt(q)*(((-1)^sgn)/A(1))*Mcd;
  else
    % Odd -- m = 1, 3, 5, 7 ...
    %fprintf('Odd Mathieu Mc(1), m = %d\n', m)
    A = mathieu_coeffs_eo(N,q,m);
    Mc = 0;
    Mcd = 0;
    for k=(N-1):-1:0
      Jks = besselj(k,s);
      Jkt = besselj(k,t);
      Jkp1s = besselj(k+1,s);
      Jkp1t = besselj(k+1,t);

      Jdks = besseljd(k,s);
      Jdkt = besseljd(k,t);
      Jdkp1s = besseljd(k+1,s);
      Jdkp1t = besseljd(k+1,t);
      
      Mc = Mc + ((-1)^k)*A(k+1)*(Jks.*Jkp1t + Jkp1s.*Jkt);
      Mcd = Mcd + ((-1)^k)*A(k+1)* ...
	    (exp(u).*(Jks.*Jdkp1t + Jkp1s.*Jdkt) - ...
	     exp(-u).*(Jdks.*Jkp1t + Jdkp1s.*Jkt) ) ;
    end
    % Do normalization.
    sgn = (m-1)/2;
    Mc = (((-1)^sgn)/A(1))*Mc;
    Mcd = sqrt(q)*(((-1)^sgn)/A(1))*Mcd;    
  end

  % At this point, fcns should be properly signed and normalized.
  
end
