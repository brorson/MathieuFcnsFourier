function plot_mathieu_coeffs_ee()
  % This plots the Fourier expansion coeffs for
  % Mathieu ce_2n.
    
  m = 50;     % Mathieu order = m-1
  % I find the peak Fourier coeff tracks m.
  N = m+10;   % I empirically find that the Fourier coeffs
              % are down by O(1e-15) about 10 positions higher
	      % than the desired order for q = 1;
  q = 100;

  % Get coeffs.
  A = mathieu_coeffs_ee(N,q,m);

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
  ce_ee = zeros(Nv,1);
  k = (0:N-1)';

  % Normalize
  %C = 2*A(1)^2 + sum(A(2:end).^2);
  %A = A/C;

  % Fourier sum
  for j=1:Nv
    ce_ee(j) = sum( A.*cos(2*k*v(j)) );
  end

  % Treat normalization of ce_0 separately
  if ((m-1) == 0)
    % Not sure why I need to scale the n=0 fcn.
    % cf.  DLMF 28.4.13.
    ce_ee = ce_ee/sqrt(2);
  end
  
  
  figure(2)
  plot(v,ce_ee)
  
end
