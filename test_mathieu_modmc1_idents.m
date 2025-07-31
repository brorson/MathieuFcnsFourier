function test_mathieu_modmc1_idents()
  % This checks modmc1 using a few identities.
    

  fail = 0;
  pass = 0;
    
  qs = [.001, .01, .1, 1, 10, 100];
  
  N = 1000;
  v = linspace(0, 10, N)';

  MM = 5;  % This is max order to test.

  %====================================================
  % Test asymptotic behavior
  fprintf('Testing asymptotic behavior per DLMF 28.20.11 ... \n')
  tol = 5e-4;
  v = linspace(3, 10, N)';  % Only want to see large v.

  % Test orders starting at m=0.
  for m=0:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      modmc1 = mathieu_modmc1(m,q,v);
      sqq = sqrt(q);
      j = besselj(m,2*sqq*cosh(v));

      %plot(v,modmc1)
      %hold on
      %plot(v,j)
      %pause()
      %close all; 
      
      % Relative norm diff
      relndiff = norm(modmc1 - j)/N;
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
  % Test Wronskian
  fprintf('Testing W(modmc1,modmc2) Wronskian per DLMF 28.20.21 ... \n')
  tol = 1e-6;
  v = linspace(0, 10, N)';
  MM = 5;  % Wronskian test starts to fail for m=6.  Must fix impls.

  % Test orders starting at m=0 for mc fcns.
  for m=0:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      [y1,y1d] = mathieu_modmc1(m,q,v);
      [y2,y2d] = mathieu_modmc2(m,q,v);

      % Compute Wronskian
      w = y1.*y2d - y1d.*y2;
      
      % Relative norm diff
      wtrue = 2/pi;
      relndiff = norm(w-wtrue)/N;
      fprintf('m = %d, q = %f, relndiff = %e ... ', m, q, relndiff)
      if (abs(relndiff) > tol)
	fprintf('Error!\n')
	fail = fail+1;
	%plot(v,w)
	%title('Wronskian')
	%pause()
	%close all; 
      else
	fprintf('\n')
	pass = pass+1;
      end
    end
  end

  fprintf('======================================\n')
  % Test Wronskian
  fprintf('Testing W(modmc1,modms2) Wronskian per DLMF 28.20.21 ... \n')
  tol = 1e-6;
  v = linspace(0, 10, N)';
  MM = 5;  % Wronskian test starts to fail for m=6.  Must fix impls.
  
  % Test orders starting at m=1.
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      [y1,y1d] = mathieu_modmc1(m,q,v);
      [y2,y2d] = mathieu_modms2(m,q,v);

      % Compute Wronskian
      w = y1.*y2d - y1d.*y2;
      
      % Relative norm diff
      wtrue = 2/pi;
      relndiff = norm(w-wtrue)/N;
      fprintf('m = %d, q = %f, relndiff = %e ... ', m, q, relndiff)
      if (abs(relndiff) > tol)
	fprintf('Error!\n')
	fail = fail+1;
	%plot(v,w)
	%title('Wronskian')
	%pause()
	%close all; 
      else
	fprintf('\n')
	pass = pass+1;
      end
    end
  end

  fprintf('At end, pass = %d, fail = %d\n', pass, fail)
  
end
