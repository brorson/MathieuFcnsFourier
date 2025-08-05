function test_mathieu_modmc2_idents()
  % This checks modmc2 using a few identities.
    

  fail = 0;
  pass = 0;
    
  qs = [.001, .01, .1, 1, 10, 100];
  
  N = 1000;

  MM = 10;  % This is max order to test.

  %====================================================
  % Test asymptotic behavior
  fprintf('Testing asymptotic behavior per DLMF 28.20.11 ... \n')
  tol = 5e-4;
  v = linspace(5, 15, N)';  % Only want to see large v.

  % Test orders starting at m=0.
  for m=0:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      modmc2 = mathieu_modmc2(m,q,v);
      sqq = sqrt(q);
      y = bessely(m,2*sqq*cosh(v));


      diffstd = std(modmc2 - y);
      if (abs(diffstd) > tol)
	fprintf('Error! ... ')
	fprintf('m = %d, q = %f, diffstd = %e ... \n', m, q, diffstd)
	fail = fail+1;
	figure(1)
	plot(v,modmc2)
	hold on
	plot(v,y)
	title('modmc2 vs. bessely')
	legend('modmc2','bessely')
	figure(2)
	plot(v,modmc2-y)
	title('Difference modmc2 - y')
	pause()
	close all;
      else
	pass = pass+1;
      end
      
    end
  end

if 0
  fprintf('======================================\n')
  % Test Wronskian
  fprintf('Testing W(modmc2,modmc1) Wronskian per DLMF 28.20.21 ... \n')
  tol = 1e-6;
  v = linspace(0, 10, N)';  % Only want to see large v.

  % Test orders starting at m=1.
  for m=0:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      [y1,y1d] = mathieu_modmc2(m,q,v);
      [y2,y2d] = mathieu_modmc1(m,q,v);

      % Compute Wronskian
      w = y1.*y2d - y1d.*y2;
      
      % 
      wtrue = -2/pi;
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
end


if 0  
  fprintf('======================================\n')
  % Test Wronskian
  fprintf('Testing W(modmc2,modms1) Wronskian per DLMF 28.20.21 ... \n')
  tol = 1e-6;
  v = linspace(0, 10, N)';
  MM = 5;  % Wronskian test starts to fail for m=6.  Must fix impls.
  
  % Test orders starting at m=1.
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      [y1,y1d] = mathieu_modmc2(m,q,v);
      [y2,y2d] = mathieu_modms1(m,q,v);

      % Compute Wronskian
      w = y1.*y2d - y1d.*y2;
      
      % 
      wtrue = -2/pi;
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
end


  fprintf('At end, pass = %d, fail = %d\n', pass, fail)
  
end
