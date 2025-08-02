function test_mathieu_modms1_idents()
  % This checks modms1 using a few identities.
    

  fail = 0;
  pass = 0;
    
  qs = [.001, .01, .1, 1, 10, 100];
  
  N = 200;
  v = linspace(0, 10, N)';

  MM = 5;  % This is max order to test.

  %====================================================
  % Test asymptotic behavior
  fprintf('Testing asymptotic behavior per DLMF 28.20.11 ... \n')
  tol = 5e-4;
  v = linspace(4, 15, N)';  % Only want to see large v.
  
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
      
      % 
      diffstd = std(modms1 - j);
      fprintf('m = %d, q = %f, diffstd = %e ... ', m, q, diffstd)
      if (abs(diffstd) > tol)
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
  fprintf('Testing W(modms1,modms2) Wronskian per DLMF 28.20.21 ... \n')
  tol = 1e-6;
  v = linspace(0, 10, N)';
  MM = 5;  % Wronskian test starts to fail for m=6.  Must fix impls.

  % Test orders starting at m=0 for mc fcns.
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      [y1,y1d] = mathieu_modms1(m,q,v);
      [y2,y2d] = mathieu_modms2(m,q,v);

      % Compute Wronskian
      w = y1.*y2d - y1d.*y2;
     
      % 
      wtrue = 2/pi;
      diffstd = std(w-wtrue);
      fprintf('m = %d, q = %f, diffstd = %e ... ', m, q, diffstd)
      if (abs(diffstd) > tol)
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
  fprintf('Testing W(modms1,modmc2) Wronskian per DLMF 28.20.21 ... \n')
  tol = 1e-6;
  v = linspace(0, 10, N)';
  MM = 5;  % Wronskian test starts to fail for m=6.  Must fix impls.
  
  % Test orders starting at m=1.
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      [y1,y1d] = mathieu_modms1(m,q,v);
      [y2,y2d] = mathieu_modmc2(m,q,v);

      % Compute Wronskian
      w = y1.*y2d - y1d.*y2;
      
      % 
      wtrue = 2/pi;
      diffstd = std(w-wtrue);
      fprintf('m = %d, q = %f, diffstd = %e ... ', m, q, diffstd)
      if (abs(diffstd) > tol)
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

  
  fprintf('At end, pass = %d, fail = %d\n', pass, fail)
  
end
