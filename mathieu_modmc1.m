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
  % This is per the book "Accurate Computation of Mathieu Functions",
  % Malcolm M. Bibby & Andrew F. Peterson.
  if ( (m>5 && q<.001) || ...
       (m>7 && q<.01) || ...
       (m>10 && q<.1) || ...
       (m>15 && q<1)  || ...
       (m>20 && q<10) || ...
       (m>30 && q<100) )
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
      if (c==0)
        % Non-adaptive calc
        Jks = besselj(k,s);
        Jkt = besselj(k,t);
        Jdks = besseljd(k,s);
        Jdkt = besseljd(k,t);
	%fprintf('Jks = %20.17e\n', Jks)
	%fprintf('Jkt = %20.17e\n', Jkt)

	% I don't use the (-1)^n idiom.  Instead I just either add or subtract
	% terms directly.  The goal is to minimize catastrophic cancellation
	% by summing pos with pos and neg with neg.  Then I subtract
	% the whole thing below.
	tt = A(k+1)*(Jks.*Jkt);
	ttd = A(k+1)*(exp(u).*Jks.*Jdkt - exp(-u).*Jdks.*Jkt);
	%fprintf('k = %d, tt = %20.17e\n', k, tt)
	
	if (mod(k,2)==0)
	  % Even terms have + sign, 
	  sgn = 1;
	else
	  % Odd terms have - sign
	  sgn = -1;
	end
	
	% Do sum using separate sums for + and -
	tt = sgn*tt;
	if (tt<0)
	  % Neg terms
	  Mcm = Mcm + tt;
	else
	  % Pos terms
	  Mcp = Mcp + tt;
	end
	
	% Do sum using separate sums for + and -
	ttd = sgn*ttd;
	if (ttd<0)
	  % Neg terms
	  Mcdm = Mcdm + ttd;
	else
	  % Pos terms
	  Mcdp = Mcdp + ttd;
	end
	
      else  % if (c==0)
        % Adaptive calc
        Jkmcs = besselj(k-c,s);
        Jkpcs = besselj(k+c,s);
        Jkpct = besselj(k+c,t);
        Jkmct = besselj(k-c,t);

        Jdkmcs = besseljd(k-c,s);
        Jdkpcs = besseljd(k+c,s);
        Jdkpct = besseljd(k+c,t);
        Jdkmct = besseljd(k-c,t);

	tt = A(k+1)*(Jkmcs.*Jkpct + Jkpcs.*Jkmct);
	ttd = A(k+1)* ...
            (exp(u).*(Jkmcs.*Jdkpct + Jkpcs.*Jdkmct) - ...
            exp(-u).*(Jdkmcs.*Jkpct + Jdkpcs.*Jkmct)) ;
	%fprintf('k = %d, tt = %e\n', k, tt)
	
	if (mod(k,2)==0)
	  % Even terms have + sign
	  sgn = 1;
	else
	  % Odd terms have - sign
	  sgn = -1;
	end

	% Do sum using separate sums for + and -
	tt = sgn*tt;
	if (tt<0)
	  Mcm = Mcm + tt;
	else
	  Mcp = Mcp + tt;
	end
	
	% Do sum using separate sums for + and -
	ttd = sgn*ttd;
	if (ttd<0)
	  Mcdm = Mcdm + ttd;
	else
	  Mcdp = Mcdp + ttd;
	end

      end
    end

    % Add pos and neg terms to get final result
    %fprintf('Mcp = %e, Mcm = %e\n', Mcp, Mcm)
    Mc = (Mcp + Mcm);
    Mcd = (Mcdp + Mcdm);

    % Do normalization.  Note normalization depends upon c.
    sgn = (-1)^(m/2);
    %fprintf('Before normalization, sgn = %d, A(c+1) = %e, Mc = %e\n', sgn, A(c+1), Mc)
    Mc = (sgn/A(c+1))*Mc;    
    Mcd = sqrt(q)*(sgn/A(c+1))*Mcd;
    %fprintf('After normalization, Mc = %e\n', Mc)

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
        % Non-adaptive calc
        Jks = besselj(k,s);
        Jkp1s = besselj(k+1,s);
        Jkt = besselj(k,t);
        Jkp1t = besselj(k+1,t);

        Jdks = besseljd(k,s);
        Jdkp1s = besseljd(k+1,s);
        Jdkt = besseljd(k,t);
        Jdkp1t = besseljd(k+1,t);
	
	%fprintf('Jks = %20.17e\n', Jks)
	%fprintf('Jkp1s = %20.17e\n', Jkp1s)
	%fprintf('Jkt = %20.17e\n', Jkt)
	%fprintf('Jkp1t = %20.17e\n', Jkp1t)	
        
	tt = A(k+1)*(Jks.*Jkp1t + Jkp1s.*Jkt);
	ttd = A(k+1)* ...
	          (exp(u).*(Jks.*Jdkp1t + Jkp1s.*Jdkt) - ...
	          exp(-u).*(Jdks.*Jkp1t + Jdkp1s.*Jkt) ) ;
	%fprintf('k = %d, tt = %20.17e\n', k, tt)
	
	if (mod(k,2)==0)
	  % Even terms have + sign
	  sgn = 1.0;
	else
	  % Odd terms have - sign
	  sgn = -1.0;
	end

	% Do sum using separate sums for + and -
	tt = sgn*tt;
	if (tt<0)
	  Mcm = Mcm + tt;
	else
	  Mcp = Mcp + tt;
	end
	
	% Do sum using separate sums for + and -
	ttd = sgn*ttd;
	if (ttd<0)
	  Mcdm = Mcdm + ttd;
	else
	  Mcdp = Mcdp + ttd;
	end

      else  % if (c==0)
        % Adaptive calc
        Jkmcs = besselj(k-c,s);
        Jkpcs = besselj(k+c+1,s);
        Jkpct = besselj(k+c+1,t);
        Jkmct = besselj(k-c,t);

        Jdkmcs = besseljd(k-c,s);
        Jdkpcs = besseljd(k+c+1,s);
        Jdkpct = besseljd(k+c+1,t);
        Jdkmct = besseljd(k-c,t);

	tt = A(k+1)*(Jkmcs.*Jkpct + Jkpcs.*Jkmct);
	ttd = A(k+1)* ...
	          (exp(u).*(Jkmcs.*Jdkpct + Jkpcs.*Jdkmct) - ...
	          exp(-u).*(Jdkmcs.*Jkpct + Jdkpcs.*Jkmct) ) ;
	%fprintf('k = %d, tt = %e, Mc = %e\n', k, tt, Mc)
	
	if (mod(k,2)==0)
	  % Even terms have + sign
	  sgn = 1.0;
	else
	  % Odd terms have - sign
	  sgn = -1.0;
	end

	% Do sum using separate sums for + and -
	tt = sgn*tt;
	if (tt<0)
	  Mcm = Mcm + tt;
	else
	  Mcp = Mcp + tt;
	end
	
	% Do sum using separate sums for + and -
	ttd = sgn*ttd;
	if (ttd<0)
	  Mcdm = Mcdm + ttd;
	else
	  Mcdp = Mcdp + ttd;
	end
	
      end

    end
    
    % Add pos and neg terms to get final result
    %fprintf('Mcp = %e, Mcm = %e\n', Mcp, Mcm)
    Mc = (Mcp + Mcm);
    %fprintf('Mcp = %20.17e, Mcm = %20.17e, Mc = %20.17e\n', Mcp, Mcm, Mc);
    Mcd = (Mcdp + Mcdm);

    % Do normalization.  Note normalization depends upon c.
    sgn = (-1)^((m-1)/2);
    %fprintf('Before normalization, sgn = %d, Mc = %e\n', sgn, Mc)
    Mc = (sgn/A(c+1))*Mc;
    Mcd = sqrt(q)*(sgn/A(c+1))*Mcd;    
    %fprintf('After normalization, Mc = %e\n', Mc)

  end

  % At this point, fcns should be properly signed and normalized.
  
end
