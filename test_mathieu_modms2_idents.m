function [pass, fail] = test_mathieu_modms2_idents()
  % This checks modms2 using a few identities.
    

  fail = 0;
  pass = 0;
    
  qs = [.001, .01, .1, 1, 10, 100];
  
  N = 1000;

  MM = 20;  % This is max order to test.

  %====================================================
  % Test asymptotic behavior
  fprintf('Testing asymptotic behavior per DLMF 28.20.11 ... \n')
  tol = 5e-4;
  v = linspace(6, 15, N)';  % Only want to see large v.

  % Test orders starting at m=1.
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      modms2 = mathieu_modms2(m,q,v);
      sqq = sqrt(q);
      y = bessely(m,2*sqq*cosh(v));

      % 
      diff = modms2-y;
      diffstd = std(diff);
      if (abs(diffstd) > tol)
	fprintf('Error! ... ')
	fprintf('m = %d, q = %f, diffstd = %e \n', m, q, diffstd)
	fail = fail+1;
	figure(1)
	plot(v,modms2)
	hold on
	plot(v,y)
	title('modms2 and bessely')
	figure(2)
	plot(v,diff)
	title('Difference')
	%pause()
	close all; 
      else
	pass = pass+1;
      end
      
    end
  end

  %====================================================  
  % Test round trip error
  NN = 100;
  v = linspace(2,5,NN)';
  h = 1e-5;
  tol = 1e-3;
  
  % Parameters to vary
  ms = 1:MM;  % Ms orders start at 1
  qs = logspace(-4,2,10);
  
  fprintf('Computing round trip error for modms2\n')
  
  % Loop over q and m
  for i=1:length(ms)
    m = ms(i);
    fprintf('-----------  m = %d  -----------\n', m)
    for j = 1:length(qs)
      q = qs(j);

      % Use sixth order deriv.  Coeffs are:
      % -1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60
      [~,fm3] = mathieu_modms2(m,q,v-3*h);
      [~,fm2] = mathieu_modms2(m,q,v-2*h);
      [~,fm1] = mathieu_modms2(m,q,v-h);
      [~,fp1] = mathieu_modms2(m,q,v+h);
      [~,fp2] = mathieu_modms2(m,q,v+2*h);
      [~,fp3] = mathieu_modms2(m,q,v+3*h);
      ydd = -fm3/60 + 3*fm2/20 - 3*fm1/4 + 3*fp1/4 - 3*fp2/20 + fp3/60;

      [y,yd] = mathieu_modms2(m,q,v);
      a = mathieu_b(m,q);

      r = ydd/(h) - (a - 2*q*cosh(2*v)).*y;
     
      stddev = std(r);
      l2norm = norm(y);
     
      if ((stddev/l2norm) > tol)
              fprintf('Error! ... ')
              fprintf('m = %d, q = %f, stddev = %e, l2norm = %e ... \n', m, q, stddev, l2norm)
              fail = fail+1;
              %figure(1)
              %plot(v,y)
              %hold on
              %plot(v,yd)
              %title('modms2 & deriv')
              %legend('y','yd')
              %figure(2)
              %plot(v,r)
              %title('Residual')
              %pause()
              %close all;   
      else
              pass = pass+1;
      end
     
    end
  end



if 1
  fprintf('======================================\n')
  % Test Wronskian
  fprintf('Testing W(modms2,modms1) Wronskian per DLMF 28.20.21 ... \n')
  tol = 1e-6;
  v = linspace(6, 15, N)';  % Only want to see large v.

  % Test orders starting at m=1.
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      [y1,y1d] = mathieu_modms2(m,q,v);
      [y2,y2d] = mathieu_modms1(m,q,v);

      % Compute Wronskian
      w = y1.*y2d - y1d.*y2;
      
      % 
      wtrue = -2/pi;
      diffstd = std(w-wtrue);
      if (abs(diffstd) > tol)
	fprintf('Error! ... ')
	fail = fail+1;
	fprintf('m = %d, q = %f, diffstd = %e\n', m, q, diffstd)
	%plot(v,w)
	%title('Wronskian')
	%pause()
	%close all; 
      else
	pass = pass+1;
      end
    end
  end
end


if 1
  fprintf('======================================\n')
  % Test Wronskian
  fprintf('Testing W(modms2,modmc1) Wronskian per DLMF 28.20.21 ... \n')
  tol = 1e-6;
  v = linspace(6, 15, N)';  % Only want to see large v.

  % Test orders starting at m=1.
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      [y1,y1d] = mathieu_modms2(m,q,v);
      [y2,y2d] = mathieu_modmc1(m,q,v);

      % Compute Wronskian
      w = y1.*y2d - y1d.*y2;
      
      % 
      wtrue = -2/pi;
      diffstd = std(w-wtrue);
      if (abs(diffstd) > tol)
	fprintf('Error! ... ')
	fail = fail+1;
	fprintf('m = %d, q = %f, diffstd = %e\n ', m, q, diffstd)
	%plot(v,w)
	%title('Wronskian')
	%pause()
	%close all; 
      else
	pass = pass+1;
      end
    end
  end
end
  
  fprintf('At end, pass = %d, fail = %d\n', pass, fail)
  
end
