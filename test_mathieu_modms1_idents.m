function test_mathieu_modms1_idents()
  % This checks modms1 using a few identities.
    

  fail = 0;
  pass = 0;
    
  qs = [.001, .01, .1, 1, 10, 100];
  
  N = 200;
  v = linspace(0, 10, N);

  MM = 35;  % This is max order to test.

  %====================================================
  % Test asymptotic behavior
  fprintf('Testing asymptotic behavior per DLMF 28.20.11 ... \n')
  tol = 5e-4;
  v = linspace(3, 10, N);  % Only want to see large v.
  
  % Test orders starting at m=1
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      modms1 = mathieu_modms1(m,q,v);
      sqq = sqrt(q);
      j = besselj(m,2*sqq*cosh(v));

      %plot(v,modms1)
      %hold on
      %plot(v,j)
      %pause()
      %close all; 
      
      % Relative norm diff
      relndiff = norm(modms1 - j)/N;
      fprintf('m = %d, q = %f, relndiff = %e ... ', m, q, relndiff)
      if (abs(relndiff) > tol)
	fprintf('Error!\n')
	fail = fail+1;
      else
	fprintf('\n')
	pass = pass+1;
      end
      
      
    end
  end

  fprintf('======================================\n')

  fprintf('At end, pass = %d, fail = %d\n', pass, fail)
  
end
