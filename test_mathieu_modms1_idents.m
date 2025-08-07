function test_mathieu_modms1_idents()
  % This checks modms1 using a few identities.
    

  fail = 0;
  pass = 0;
    
  qs = [.001, .01, .1, 1, 10, 100];
  
  N = 200;  % Number of sample pts in v domain.
  MM = 20;  % This is max order to test.

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

      diffstd = std(modms1 - j);

      if (abs(diffstd) > tol)
	fprintf('Error! ... ')
	fprintf('m = %d, q = %f, diffstd = %e ... \n', m, q, diffstd)
	fail = fail+1;
	%figure(1)
	%plot(v,modms1)
	%hold on
	%plot(v,j)
	%title('modms1 vs. besselj')
	%legend('modms1','besselj')
	%figure(2)
	%plot(v,modms1-j)
	%title('Difference modms1 - j')
	%pause()
	%close all;  
      else
	pass = pass+1;
      end
     
    end
  end


  %====================================================  
  % Test round trip error
  NN = 100;
  v = linspace(2,5,NN)';
  h = 1e-4;
  tol = 1e-3;
  
  % Parameters to vary
  ms = 1:MM;  % Ms orders start at 1
  qs = logspace(-4,2,10);
  
  fprintf('Computing round trip error for modms1\n')
  
  % Loop over q and m
  for i=1:length(ms)
    m = ms(i);
    fprintf('-----------  m = %d  -----------\n', m)
    for j = 1:length(qs)
      q = qs(j);

      % Use sixth order deriv.  Coeffs are:
      % -1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60
      [~,fm3] = mathieu_modms1(m,q,v-3*h);
      [~,fm2] = mathieu_modms1(m,q,v-2*h);
      [~,fm1] = mathieu_modms1(m,q,v-h);
      [~,fp1] = mathieu_modms1(m,q,v+h);
      [~,fp2] = mathieu_modms1(m,q,v+2*h);
      [~,fp3] = mathieu_modms1(m,q,v+3*h);
      ydd = -fm3/60 + 3*fm2/20 - 3*fm1/4 + 3*fp1/4 - 3*fp2/20 + fp3/60;

      [y,yd] = mathieu_modms1(m,q,v);
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
	      %title('modms1 & deriv')
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

 
if 0
  %====================================================  
  fprintf('======================================\n')
  % Test Wronskian
  fprintf('Testing W(modms1,modms2) Wronskian per DLMF 28.20.21 ... \n')
  tol = 1e-6;
  v = linspace(0, 10, N)';
  MM = 5;  % Wronskian test starts to fail for m=6.  Must fix impls.

  % Test orders starting at m=1 for ms fcns.
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
end


if 0
  %==================================================== 
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
end

  
  fprintf('At end, pass = %d, fail = %d\n', pass, fail)
  
end
