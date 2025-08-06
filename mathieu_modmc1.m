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
 
  % Set offset used in Bessel fcn depending upon order m.
  % This is per the book "Accurately Calculating Mathieu Functions",
  % XXXXX & YYYY
  if ( (m>10 && q<.1) || (m>20 && q<100) )
    c = floor(m/2);
  else
    c = 0;
  end
%c = 0;

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
    % Even -- m = 0, 2, 4, 6, ...
    %fprintf('Even Mathieu Mc(1), m = %d\n', m)
    A = mathieu_coeffs_ee(N,q,m);
    Mcp = zeros(size(u));
    Mcm = zeros(size(u));
    Mcdp = zeros(size(u));
    Mcdm = zeros(size(u));
    for k=(N-1):-1:0
    %for k=0:(N-1)
      if (c==0)
        Jks = besselj(k,s);
        Jkt = besselj(k,t);

        Jdks = besseljd(k,s);
        Jdkt = besseljd(k,t);

	% I don't use the (-1)^n idiom.  Instead I just either add or subtract
	% terms directly.  The goal is to minimize catastrophic cancellation
	% by summing pos with pos and neg with neg.  Then I subtract
	% the whole thing below.
        if (mod(k,2) == 0)
          % Additive terms
          Mcp = Mcp + A(k+1)*(Jks.*Jkt);
          Mcdp = Mcdp + A(k+1)*(exp(u).*Jks.*Jdkt - exp(-u).*Jdks.*Jkt);
         else
          % Subtractive terms.
          Mcm = Mcm + A(k+1)*(Jks.*Jkt);
          Mcdm = Mcdm + A(k+1)*(exp(u).*Jks.*Jdkt - exp(-u).*Jdks.*Jkt);
        end

      else
        Jkmcs = besselj(k-c,s);
        Jkpcs = besselj(k+c,s);
        Jkpct = besselj(k+c,t);
        Jkmct = besselj(k-c,t);

        Jdkmcs = besseljd(k-c,s);
        Jdkpcs = besseljd(k+c,s);
        Jdkpct = besseljd(k+c,t);
        Jdkmct = besseljd(k-c,t);

        if (mod(k,2) == 0)
          % Even plus terms
          Mcp = Mcp + A(k+1)*(Jkmcs.*Jkpct + Jkpcs.*Jkmct);
          Mcdp = Mcdp + A(k+1)* ...
            (exp(u).*(Jkmcs.*Jdkpct + Jkpcs.*Jdkmct) - ...
            exp(-u).*(Jdkmcs.*Jkpct + Jdkpcs.*Jkmct)) ;
        else
          % Odd negative terms
          Mcm = Mcm + A(k+1)*(Jkmcs.*Jkpct + Jkpcs.*Jkmct);
          Mcdm = Mcdm + A(k+1)* ...
            (exp(u).*(Jkmcs.*Jdkpct + Jkpcs.*Jdkmct) - ...
            exp(-u).*(Jdkmcs.*Jkpct + Jdkpcs.*Jkmct)) ;
        end
      end
    end
    % Subtract for now.
    % Later implement sort and sum, or Kahan summation if needed.
    Mc = Mcp - Mcm;
    %fprintf('Mcp = %e, Mcm = %e, Mc = %e\n', Mcp, Mcm, Mc)    
    Mcd = Mcdp - Mcdm;

    % Do normalization.  Note normalization depends upon c.
    sgn = m/2;
    Mc = (((-1)^sgn)/A(c+1))*Mc;    
    Mcd = sqrt(q)*(((-1)^sgn)/A(c+1))*Mcd;

  else
    % Odd -- m = 1, 3, 5, 7 ...
    %fprintf('Odd Mathieu Mc(1), m = %d\n', m)
    A = mathieu_coeffs_eo(N,q,m);
    Mcp = zeros(size(u));
    Mcm = zeros(size(u));
    Mcdp = zeros(size(u));
    Mcdm = zeros(size(u));
    for k=(N-1):-1:0
      if (c==0)
        Jks = besselj(k,s);
        Jkp1s = besselj(k+1,s);
        Jkt = besselj(k,t);
        Jkp1t = besselj(k+1,t);
  
        Jdks = besseljd(k,s);
        Jdkp1s = besseljd(k+1,s);
        Jdkt = besseljd(k,t);
        Jdkp1t = besseljd(k+1,t);
        
        if (mod(k,2) == 0)
          % Even plus terms
          Mcp = Mcp + A(k+1)*(Jks.*Jkp1t + Jkp1s.*Jkt);
          Mcdp = Mcdp + A(k+1)* ...
	          (exp(u).*(Jks.*Jdkp1t + Jkp1s.*Jdkt) - ...
	          exp(-u).*(Jdks.*Jkp1t + Jdkp1s.*Jkt) ) ;
        else
          % Odd minus terms.
          Mcm = Mcm + A(k+1)*(Jks.*Jkp1t + Jkp1s.*Jkt);
          Mcdm = Mcdm + A(k+1)* ...
	          (exp(u).*(Jks.*Jdkp1t + Jkp1s.*Jdkt) - ...
	          exp(-u).*(Jdks.*Jkp1t + Jdkp1s.*Jkt) ) ;
        end

      else
        Jkmcs = besselj(k-c,s);
        Jkpcs = besselj(k+c+1,s);
        Jkpct = besselj(k+c+1,t);
        Jkmct = besselj(k-c,t);

        Jdkmcs = besseljd(k-c,s);
        Jdkpcs = besseljd(k+c+1,s);
        Jdkpct = besseljd(k+c+1,t);
        Jdkmct = besseljd(k-c,t);

        if (mod(k,2) == 0)
          % Even plus terms
          Mcp = Mcp + A(k+1)*(Jkmcs.*Jkpct + Jkpcs.*Jkmct);
          Mcdp = Mcdp + A(k+1)* ...
	          (exp(u).*(Jkmcs.*Jdkpct + Jkpcs.*Jdkmct) - ...
	          exp(-u).*(Jdkmcs.*Jkpct + Jdkpcs.*Jkmct) ) ;
        else
          % Odd minus terms
          Mcm = Mcm + A(k+1)*(Jkmcs.*Jkpct + Jkpcs.*Jkmct);
          Mcdm = Mcdm + A(k+1)* ...
	          (exp(u).*(Jkmcs.*Jdkpct + Jkpcs.*Jdkmct) - ...
	          exp(-u).*(Jdkmcs.*Jkpct + Jdkpcs.*Jkmct) ) ;
        end
      end

    end
    Mc = Mcp - Mcm;
    %fprintf('Mcp = %e, Mcm = %e, Mc = %e\n', Mcp, Mcm, Mc)    
    Mcd = Mcdp - Mcdm;

    % Do normalization.  Note normalization depends upon c.
    sgn = (m-1)/2;
    Mc = (((-1)^sgn)/A(c+1))*Mc;
    Mcd = sqrt(q)*(((-1)^sgn)/A(c+1))*Mcd;    
  end

  % At this point, fcns should be properly signed and normalized.
  
end
