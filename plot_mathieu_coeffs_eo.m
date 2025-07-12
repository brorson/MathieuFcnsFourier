function plot_mathieu_coeffs_eo()
  % This plots the Fourier expansion coeffs for
  % Mathieu ce_2n+1.
    
  m = 3;     % Mathieu order = m-1
  % I find the peak Fourier coeff tracks m.
  N = m+10;   % I empirically find that the Fourier coeffs
              % are down by O(1e-15) about 10 positions higher
	      % than the desired order for q = 1;
  q = 1;

  % Get coeffs.
  A = mathieu_coeffs_eo(N,q,m);

  % By defintion, all fcns are positive for v = 0.  Modify coeffs
  % to obey this definition.
  if (A(1)<0)
    A = -A;
  end
  
  figure(1)
  nn = 1:length(A);
  semilogy(nn,abs(A),'bo');

  % Do check of result by plotting the fcn itself.
  Nv = 100;
  v = linspace(0,pi/2,Nv);
  ce_eo = zeros(Nv,1);
  k = (0:N-1)';

  % Normalize
  %C = 2*A(1)^2 + sum(A(2:end).^2);
  %A = A/C;

  % Fourier sum
  for j=1:Nv
    ce_eo(j) = sum( A.*cos((2*k+1)*v(j)) );
  end

  figure(2)
  plot(v,ce_eo)
  
end
