function [Mc,Mcd] = mathieu_modmc2(m, q, u)
  % This computes the modified Mathieu fcn of the 
  % second kind of order m for frequency
  % parameter q and radius u.  This is the fcn frequently
  % called Mc in the literature.
  % Inputs m, q must be scalars, u may be a vector.
    
  % Here I use the expressions given in Zhang and Jin.


  if (q<0)
    error('Modified Mathieu fcns for negatative q not implemented yet!\n')
  end

  % Force v to be col vector so returns are col vectors.
  if (size(u, 2)>1)
    u = u';
  end

  % Set offset used in Bessel fcn depending upon order m.
  % This is per the book  "Accurate Computation of Mathieu Functions",
  % Malcolm M. Bibby & Andrew F. Peterson.
  if ( (m>10 && q<.1) || (m>20 && q<100) )
    c = floor(m/2);
  else
    c = 0;
  end
%  c = 0;

  % I find the peak Fourier coeff tracks m.  Therefore
  % I adjust the matrix size based on order m.
  N = m+25;

  % Utility vars.
  sqq = sqrt(q);
  s = sqq*exp(-u);
  t = sqq*exp(u);

  % Use different coeffs depending upon whether m is even or
  % odd.
  tol = 1e-14;
  if (abs(mod(m,2)) < tol)
    % Even
    %fprintf('Even Mathieu Mc(2), m = %d\n', m)
    A = mathieu_coeffs_ee(N,q,m);
    Mcp = zeros(size(u));
    Mcm = zeros(size(u));
    Mcdp = zeros(size(u));
    Mcdm = zeros(size(u));
    for k=(N-1):-1:0
      if (c==0)
        Jks = besselj(k,s);
        Ykt = bessely(k,t);

        Jdks = besseljd(k,s);
        Ydkt = besselyd(k,t);
        
        if (mod(k,2) == 0)
          % Even plus terms
          Mcp = Mcp + A(k+1)*Jks.*Ykt ;
          Mcdp = Mcdp + A(k+1)*(exp(u).*Jks.*Ydkt - exp(-u).*Jdks.*Ykt) ;
        else
          % Odd minus terms
          Mcm = Mcm + A(k+1)*Jks.*Ykt ;
          Mcdm = Mcdm + A(k+1)*(exp(u).*Jks.*Ydkt - exp(-u).*Jdks.*Ykt) ;
        end
      else
        Jkmcs = besselj(k-c,s);
        Jkpcs = besselj(k+c,s);
        Ykpct = bessely(k+c,t);
        Ykmct = bessely(k-c,t);

        Jdkmcs = besseljd(k-c,s);
        Jdkpcs = besseljd(k+c,s);
        Ydkpct = besselyd(k+c,t);
        Ydkmct = besselyd(k-c,t);
        if (mod(k,2) == 0)
          % Even plus terms
          Mcp = Mcp + A(k+1)*(Jkmcs.*Ykpct + Jkpcs.*Ykmct);
          Mcdp = Mcdp + A(k+1)* ...
            (exp(u).*(Jkmcs.*Ydkpct + Jkpcs.*Ydkmct) - ...
             exp(-u).*(Jdkmcs.*Ykpct + Jdkpcs.*Ykmct)) ;
        else
          % Odd minus terms
          Mcm = Mcm + A(k+1)*(Jkmcs.*Ykpct + Jkpcs.*Ykmct);
          Mcdm = Mcdm + A(k+1)* ...
            (exp(u).*(Jkmcs.*Ydkpct + Jkpcs.*Ydkmct) - ...
             exp(-u).*(Jdkmcs.*Ykpct + Jdkpcs.*Ykmct)) ;
        end

      end
    end

    % Later do sort and sum, or Kahan summation
    Mc = Mcp - Mcm;
    %fprintf('Mcp = %e, Mcm = %e, Mc = %e\n', Mcp, Mcm, Mc) 
    Mcd = Mcdp - Mcdm;

    % Do normalization.  Note normalization depends upon c.
    sgn = m/2;
    Mc = (((-1)^sgn)/A(c+1))*Mc;
    Mcd = sqrt(q)*(((-1)^sgn)/A(c+1))*Mcd;

  else
    % Odd
    %fprintf('Odd Mathieu Mc(2), m = %d\n', m)
    A = mathieu_coeffs_eo(N,q,m);
    Mcp = zeros(size(u));
    Mcm = zeros(size(u));
    Mcdp = zeros(size(u));
    Mcdm = zeros(size(u));
    for k=(N-1):-1:0
      if (c==0)
        Jks = besselj(k,s);
        Ykt = bessely(k,t);
        Jkp1s = besselj(k+1,s);
        Ykp1t = bessely(k+1,t);
  
        Jdks = besseljd(k,s);
        Ydkt = besselyd(k,t);
        Jdkp1s = besseljd(k+1,s);
        Ydkp1t = besselyd(k+1,t);

        if (mod(k,2) == 0)
          % Even plus terms
          Mcp = Mcp + A(k+1)*(Jks.*Ykp1t + Jkp1s.*Ykt);
          Mcdp = Mcdp + A(k+1)* ...
              (exp(u).*(Jks.*Ydkp1t + Jkp1s.*Ydkt) - ...
               exp(-u).*(Jdks.*Ykp1t + Jdkp1s.*Ykt) ) ;
        else
          % Odd minus terms
          Mcm = Mcm + A(k+1)*(Jks.*Ykp1t + Jkp1s.*Ykt);
          Mcdm = Mcdm + A(k+1)* ...
              (exp(u).*(Jks.*Ydkp1t + Jkp1s.*Ydkt) - ...
               exp(-u).*(Jdks.*Ykp1t + Jdkp1s.*Ykt) ) ;
        end

      else
        Jkmcs = besselj(k-c,s);
        Jkpcs = besselj(k+c+1,s);
        Ykpct = bessely(k+c+1,t);
        Ykmct = bessely(k-c,t);
  
        Jdkmcs = besseljd(k-c,s);
        Jdkpcs = besseljd(k+c+1,s);
        Ydkpct = besselyd(k+c+1,t);
        Ydkmct = besselyd(k-c,t);

        if (mod(k,2) == 0)
          % Even plus terms
          Mcp = Mcp + A(k+1)*(Jkmcs.*Ykpct + Jkpcs.*Ykmct);
          Mcdp = Mcdp + A(k+1)* ...
                (exp(u).*(Jkmcs.*Ydkpct + Jkpcs.*Ydkmct) - ...
                 exp(-u).*(Jdkmcs.*Ykpct + Jdkpcs.*Ykmct) ) ;
        else
          % Odd minus terms
          Mcm = Mcm + A(k+1)*(Jkmcs.*Ykpct + Jkpcs.*Ykmct);
          Mcdm = Mcdm + A(k+1)* ...
                (exp(u).*(Jkmcs.*Ydkpct + Jkpcs.*Ydkmct) - ...
                 exp(-u).*(Jdkmcs.*Ykpct + Jdkpcs.*Ykmct) ) ;
        end

      end
    end

    % Later do sort and sum, or Kahan summation
    Mc = Mcp - Mcm;
    %fprintf('Mcp = %e, Mcm = %e, Mc = %e\n', Mcp, Mcm, Mc) 
    Mcd = Mcdp - Mcdm;

    % Do normalization.  Note normalization depends upon c.
    sgn = (m-1)/2;
    Mc = (((-1)^sgn)/A(c+1))*Mc;    
    Mcd = sqrt(q)*(((-1)^sgn)/A(c+1))*Mcd;
  end

  % At this point, fcns should be properly normalized.
  
end
