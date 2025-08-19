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
  % This is per the book "Accurate Computation of Mathieu Functions",
  % Malcolm M. Bibby & Andrew F. Peterson.
%  if ( (m>10 && q<.1) || (m>20 && q<100) )
  if ( (m>5 && q<.001) || ...
       (m>7 && q<.01) || ...
       (m>10 && q<.1) || ...
       (m>15 && q<1)  || ...
       (m>20 && q<10) || ...
       (m>30 && q<100) )
    c = floor((m-1)/2);
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
  if (abs(mod(m,2) < tol))
    % Even
    %fprintf('Even Mathieu Ms(2), m = %d\n', m)
    B = mathieu_coeffs_oe(N,q,m);
    Msp = zeros(size(u));
    Msm = zeros(size(u));
    Msdp = zeros(size(u));
    Msdm = zeros(size(u));

    for k=(N-1):-1:0
      if (c==0)
	% Non-adaptive calc
	Jks = besselj(k,s);
	Ykp2t = bessely(k+2,t);
	Jkp2s = besselj(k+2,s);
	Ykt = bessely(k,t);

	Jdks = besseljd(k,s);
	Ydkp2t = besselyd(k+2,t);
	Jdkp2s = besseljd(k+2,s);
	Ydkt = besselyd(k,t);

	tt =  B(k+1).*(Jks.*Ykp2t - Jkp2s.*Ykt);
	ttd =  B(k+1)*...
		 (exp(u).*(Jks.*Ydkp2t - Jkp2s.*Ydkt) - ...
		  exp(-u).*(Jdks.*Ykp2t - Jdkp2s.*Ykt));
	
        if (mod(k,2) == 0)
	  sgn = 1;
	else
	  sgn = -1;
	end

	% Do sum using separate sums for + and -
	tt = sgn*tt;
	%Msp = Msp + tt;
	if (tt<0) 
	  % Neg terms
	  Msm = Msm + tt;
	else 
	  % Pos terms
	  Msp = Msp + tt;
	end
      
	% Do sum using separate sums for + and -
	ttd = sgn*ttd;
	%Msdp = Msdp + tt;	
	if (ttd<0) 
	  % Neg terms
	  Msdm = Msdm + ttd;
	else 
	  % Pos terms
	  Msdp = Msdp + ttd;
	end
      
      else  % if (c==0) ...
	% Adaptive calc.
	Jkmcs = besselj(k-c,s);
        Ykpct = bessely(k+c+2,t);
        Jkpcs = besselj(k+c+2,s);
        Ykmct = bessely(k-c,t);
	
        Jdkmcs = besseljd(k-c,s);
        Ydkpct = besselyd(k+c+2,t);
        Jdkpcs = besseljd(k+c+2,s);
        Ydkmct = besselyd(k-c,t);
	
	tt = B(k+1)*(Jkmcs.*Ykpct - Jkpcs.*Ykmct);
	ttd =  B(k+1)* ...
		 (exp(u).*(Jkmcs.*Ydkpct - Jkpcs.*Ydkmct) - ...
		  exp(-u).*(Jdkmcs.*Ykpct - Jdkpcs.*Ykmct));
	
        if (mod(k,2) == 0)
	  sgn = 1;
	else
	  sgn = -1;
	end

	% Do sum using separate sums for + and -
	tt = sgn*tt;
	%Msp = Msp + tt;
	if (tt<0) 
	  % Neg terms
	  Msm = Msm + tt;
	else 
	  % Pos terms
	  Msp = Msp + tt;
	end
      
	% Do sum using separate sums for + and -
	ttd = sgn*ttd;
	%Msdp = Msdp + tt;
	if (ttd<0) 
	  % Neg terms
	  Msdm = Msdm + ttd;
	else 
	  % Pos terms
	  Msdp = Msdp + ttd;
	end

      end  % if (c==0) ...
      
    end   % for k = ...
    % Sum pos and neg terms to get final answer
    Ms = Msp + Msm;
    %fprintf('Msp = %e, Msm = %e, Ms = %e\n', Msp, Msm, Ms)    
    Msd = Msdp + Msdm;
    
    % Do normalization.
    sgn = (m-2)/2;
    Ms = (((-1)^sgn)/B(c+1))*Ms;    
    Msd = sqrt(q)*(((-1)^sgn)/B(c+1))*Msd; 
    
    
  else
    % Odd
    %fprintf('Odd Mathieu Ms(2), m = %d\n', m)
    B = mathieu_coeffs_oo(N,q,m);
    Msp = zeros(size(u));
    Msm = zeros(size(u));
    Msdp = zeros(size(u));
    Msdm = zeros(size(u));  
    
    for k=(N-1):-1:0   
      if (c==0)
	% Non-adaptive
	Jks = besselj(k,s);
	Ykt = bessely(k,t);
	Jkp1s = besselj(k+1,s);
	Ykp1t = bessely(k+1,t);
	
	Jdks = besseljd(k,s);
	Ydkt = besselyd(k,t);
	Jdkp1s = besseljd(k+1,s);
	Ydkp1t = besselyd(k+1,t);
	
	tt = B(k+1).*(Jks.*Ykp1t - Jkp1s.*Ykt);
	ttd =  B(k+1)*...
		 (exp(u).*(Jks.*Ydkp1t - Jkp1s.*Ydkt) - ...
		  exp(-u).*(Jdks.*Ykp1t - Jdkp1s.*Ykt)); 

        if (mod(k,2) == 0)
	  sgn = 1;
	else
	  sgn = -1;
	end

	% Do sum using separate sums for + and -
	tt = sgn*tt;
	%Msp = Msp + tt;
	if (tt<0) 
	  % Neg terms
	  Msm = Msm + tt;
	else 
	  % Pos terms
	  Msp = Msp + tt;
	end
      
	% Do sum using separate sums for + and -
	ttd = sgn*ttd;
	%Msdp = Msdp + tt;
	if (ttd<0) 
	  % Neg terms
	  Msdm = Msdm + ttd;
	else 
	  % Pos terms
	  Msdp = Msdp + ttd;
	end

	
      else % if (c==0)
	% Adaptive calc
	Jkmcs = besselj(k-c,s);
        Jkpcs = besselj(k+c+1,s);
        Ykpct = bessely(k+c+1,t);
        Ykmct = bessely(k-c,t);
	
        Jdkmcs = besseljd(k-c,s);
        Jdkpcs = besseljd(k+c+1,s);
        Ydkpct = besselyd(k+c+1,t);
        Ydkmct = besselyd(k-c,t);
	
	tt = B(k+1)*(Jkmcs.*Ykpct - Jkpcs.*Ykmct);
	ttd = B(k+1)* ...
		 (exp(u).*(Jkmcs.*Ydkpct - Jkpcs.*Ydkmct) - ...
                  exp(-u).*(Jdkmcs.*Ykpct - Jdkpcs.*Ykmct) ) ; 

        if (mod(k,2) == 0)
	  sgn = 1;
	else
	  sgn = -1;
	end

	% Do sum using separate sums for + and -
	tt = sgn*tt;
	%Msp = Msp + tt;
	if (tt<0) 
	  % Neg terms
	  Msm = Msm + tt;
	else 
	  % Pos terms
	  Msp = Msp + tt;
	end
      
	% Do sum using separate sums for + and -
	ttd = sgn*ttd;
	%Msdp = Msdp + tt;
	if (ttd<0) 
	  % Neg terms
	  Msdm = Msdm + ttd;
	else 
	  % Pos terms
	  Msdp = Msdp + ttd;
	end

      end

    end  % for k = ...

    % Sum pos and neg terms to get final answer
    Ms = Msp + Msm;
    %fprintf('Msp = %e, Msm = %e, Ms = %e\n', Msp, Msm, Ms)    
    Msd = Msdp + Msdm;

    % Do normalization.
    sgn = (m-1)/2;
    Ms = (((-1)^sgn)/B(c+1))*Ms;
    Msd = sqrt(q)*(((-1)^sgn)/B(c+1))*Msd;     
  end

  % At this point, fcns should be properly normalized.
  
end
