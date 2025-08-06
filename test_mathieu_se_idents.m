function test_mathieu_se_idents()
  % This checks se using a few identities.

  pass = 0;
  fail = 0;
    
  qs = logspace(-3,3,21);
  qs = [-qs, qs];
  %qs = [-1000, -100, -10, -1, -0.1, -0.01, -0.001, 0, 0.0001, .001, .01, .1, 1, 10, 100, 1000];
  %qs = [1, 10];    
  
  N = 1000;
  v = linspace(-pi,pi,N)';

  MM = 10;   % This is max order to test -- usually 35
  
  
  %====================================================
  % First test normalization per DLMF 28.2.30
  fprintf('Testing normalization DLMF 28.2.30 ... \n')
  MM = 25;   % This is max order to test
  tol = 1e-14;
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      se = mathieu_se(m,q,v);
      s = trapz(v,se.*se);
      
      diff = s - pi;
      %fprintf('s = %f, err = %e\n', s, diff)
      if (abs(diff) > tol)
	      fprintf('Error!  m = %d, q = %f, diff = %e\n', m, q, diff)
	      fail = fail+1;
      else
	pass = pass+1;
      end
      
    end
  end

  fprintf('======================================\n')

  
  %====================================================
if 0
  % Next test orthogonality per DLMF 28.2,32
  fprintf('Testing orthogonality DLMF 28.2,32 ... \n')
  tol = 1e-11;
  %MM = 10;  % Max order to test
  for m1=1:MM;   for m2=m1:MM
    if (m1 == m2)
      continue
    end
    fprintf('-----------  [m1,m2] = [%d,%d]  -----------\n', m1,m2)
    for i = 1:length(qs)
      q = qs(i);
      
      % integral is more recommended over trapz
      % For some reason, I need to transpose the output of mathieu_ce.      
      se = @(m,q,v) mathieu_se(m,q,v)'; 
      f = @(v) se(m1,q,v).*se(m2,q,v);
      s = integral(f, -pi, pi);
      
      diff = s;
      %fprintf('s = %f, err = %e\n', s, diff)
      if (abs(diff) > tol)
	      fprintf('Error!  [m1,m2] = [%d,%d], q = %f, diff = %e\n', m1, m2, q, diff)
	      fail = fail+1;
      else
	pass = pass+1;
      end
    end
  end; end
  
  fprintf('======================================\n')
end

  %====================================================
if 1
  % Test high order finite diff
  fprintf('Testing finite diff ... \n')
  tol = 1e-6;
  NN = 5000;
  v = linspace(-pi,pi*(NN-1)/NN,NN)';
  MM = 10;
  qs = [-100, -10, -1, -0.1, -0.01, -0.001, 0, 0.0001, .001, .01, .1, 1, 10, 100];

  % Test orders starting at m=0 for mc fcns.
  h = v(2)-v(1);
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    for i = 1:length(qs)
      q = qs(i);
      
      [y,yd] = mathieu_se(m,q,v);
      a = mathieu_b(m,q);

      ydd = fd_deriv(yd);
      r = ydd/(h) + (a - 2*q*cos(2*v)).*y;

    % Err -- fix this later.
      diffstd = std(r(5:(end-4)));
      %fprintf('m = %d, q = %f, diffstd = %e ... ', m, q, diffstd)
      if (abs(diffstd) > tol)
        fprintf('Error! ... ')
        fprintf('m = %d, q = %f, diffstd = %e\n', m, q, diffstd)
        fail = fail+1;
        %figure(1)
        %plot(v,y);
        %hold on
        %plot(v,yd);
        %title('fcn')
        %legend('y','yd')
        %figure(3)
        %plot(v(3:(end-2)),r(3:(end-2)))
        %title('Round trip residual')
        %pause()
        %close all; 
      else
        %fprintf('\n')
        pass = pass+1;
      end
    end
  end
  
  fprintf('======================================\n')

end

  
  %====================================================
  % Test q = 0 case per DLMF 28.2.29
  fprintf('Test se tends to sin for q = -1e-13 ... \n')

  N = 1000;
  v = linspace(-pi,pi,N)';

  tol = 2e-13;
  q = -1e-13;
  MM = 10;  % Max order to test
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    LHS = mathieu_se(m,q,v);
    RHS = sin(m*v);
    
    diffstd = std(LHS-RHS);
    if (diffstd > tol)
      fprintf('Error!  m = %d, q = %5.3f, LHS(floor(N/2)) = %f, RHS(floor(N/2)) = %f, diffstd = %e\n', ...
	      m, q, LHS(1), RHS(1), diffstd)
      fail = fail+1;
      figure(1)
      plot(v,LHS);
      hold on
      plot(v,RHS);
      title('fcn')
      legend('se','sin')
      pause()
      close all;
    else
      pass = pass+1;
    end
  end
  
  fprintf('Test se tends to sin for q = 0 ... \n')
  q = 0;
  %MM = 10;  % Max order to test
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    LHS = mathieu_se(m,q,v);
    RHS = sin(m*v);
    
    diffstd = std(LHS-RHS);
    %fprintf('m = %d, LHS(1) = %f, RHS(1) = %f, std = %e\n', m, LHS(1), RHS(1), diffstd)
    if (diffstd > tol)
      fprintf('Error!  m = %d, q = %5.3f, LHS(1) = %f, RHS(1) = %f, diffstd = %e\n', m, q, LHS(1), RHS(1), diffstd)
      fail = fail+1;
    else
      pass = pass+1;
    end
  end
  
    fprintf('Test se tends to sin for q = 1e-13 ... \n')
  q = 1e-13;
  %MM = 10;  % Max order to test
  for m=1:MM
    fprintf('-----------  m = %d  -----------\n', m)
    LHS = mathieu_se(m,q,v);
    RHS = sin(m*v);
    
    diffstd = std(LHS-RHS);
    %fprintf('m = %d, LHS(1) = %f, RHS(1) = %f, std = %e\n', m, LHS(1), RHS(1), diffstd)
    if (diffstd > tol)
      fprintf('Error!  m = %d, q = %5.3f, LHS(1) = %f, RHS(1) = %f, diffstd = %e\n', m, q, LHS(1), RHS(1), diffstd)
      fail = fail+1;
    else
      pass = pass+1;
    end
  end
  
  fprintf('======================================\n')
  
  %====================================================
  % Next test even fcns per DLMF 28.2,34
  fprintf('Test rotation identity for even fcns per DLMF 28.2,34 ... \n')
  tol = 1e-13;
  v = 0.05;  % Just check one random point.
  %MM = 10;  % Sets max order to test
  for m=2:2:MM
    ss = ((-1)^((m-2)/2));
    fprintf('-----------  m = %d, ss = %d  -----------\n', m, ss)
    %fprintf('s = %f\n', s)
    for i = 1:length(qs)
      q = qs(i);
      LHS = mathieu_se(m,-q,v);
      RHS = ss*mathieu_se(m,q,(pi/2)-v);
      diff = LHS-RHS;
      %fprintf('m = %d, LHS = %f, RHS = %f, err = %e\n', m, LHS, RHS, diff)
      if (abs(diff) > tol)
	      fprintf('Error!  m = %d, q = %5.3f, LHS = %f, RHS = %f, diff = %e\n', m, q, LHS, RHS, diff)
	      fail = fail+1;
      else
	pass = pass+1;
      end
    end
  end
  
  fprintf('======================================\n')  
 
  %====================================================
  % Next test odd fcns per DLMF 28.2,35
  fprintf('Test rotation identity for odd fcns per DLMF 28.2,35 ... \n')
  v = 0.05;  % Just check one random point.
  %MM = 10;  % Sets max order to test
  for m=1:2:MM+1
    ss = ((-1)^((m-1)/2));
    fprintf('-----------  m = %d, ss = %d  -----------\n', m, ss)
    %fprintf('s = %f\n', s)
    for i = 1:length(qs)
      q = qs(i);
      LHS = mathieu_se(m,-q,v);
      RHS = ss*mathieu_ce(m,q,pi/2-v);
      diff = LHS-RHS;
      %fprintf('m = %d, LHS = %f, RHS = %f, err = %e\n', m, LHS, RHS, diff)
      if (abs(diff) > tol)
	      fprintf('Error!  m = %d, q = %5.3f, LHS = %f, RHS = %f, diff = %e\n', m, q, LHS, RHS, diff)
	      fail = fail+1;
      else
	pass = pass+1;
      end
    end
  end
  fprintf('======================================\n')  
 
if 0
  %====================================================
  % Next test small q expansions per DLMF 28.6.23.  The goal
  % is to make sure my fcn impls go the right way for negative
  % q values.
  fprintf('Test small q expansions per DLMF 28.6.23 ... \n')

  tol = 1e-4;
  N = 1000;
  v = linspace(-pi,pi,N)';
  
  qs = [-1.5, -1, -.5, -.2, -.1, .1, .2, .5, 1, 1.5];

  MM = 35;

  % The first three orders don't work unless q is really small
  %-----------------------------------------------------------
  % m = 1
  qs = [-.2, -.1, .1, .2];  
  for i=1:length(qs)
    q = qs(i);
    % Define power series after setting q.
    se1 = @(v) sin(v) - q*sin(3*v)/8 + (q^2)*(2*sin(5*v)/3 + 2*sin(3*v) - sin(v))/128 ...
	  - (q^3)*(sin(7*v)/9 + 8*sin(5*v)/9 - sin(3*v)/3 - 2*sin(v))/1024;
    m = 1;
    fprintf('----------- m = %d, q = %f  -----------\n', m, q)  
    dlmfse1 = se1(v);
    myse1 = mathieu_se(m,q,v);
    %plot(v, dlmfse1,'b-')
    %hold on
    %plot(v, myse1,'r.')
    %legend('dlmf','me')
    diffstd = std(dlmfse1 - myse1);
    if (diffstd > tol)
      fprintf('Error!  m = %d, q = %5.3f, diffstd = %e\n', m, q, diffstd)
      fail = fail+1;
    else
      pass = pass+1;
    end
    %pause()
    %close all;
  end

  %-----------------------------------------------------------  
  % m = 2
  qs = [-.5, -.2, -.1, .1, .2, .5];    
  for i=1:length(qs)
    % Define power series after setting q.
    se2 = @(v) sin(2*v) - q*sin(4*v)/12 + (q^2)*(sin(6*v)/3 - 4*sin(2*v)/9)/128;
    m = 2;
    fprintf('----------- m = %d, q = %f  -----------\n', m, q)  
    dlmfse2 = se2(v);
    myse2 = mathieu_se(m,q,v);
    %plot(v, dlmfse2)
    %hold on
    %plot(v, myse2)
    %legend('dlmf','me')
    diffstd = std(dlmfse2 - myse2);
    if (diffstd > tol)
      fprintf('Error!  m = %d, q = %5.3f, diffstd = %e\n', m, q, diffstd)
      fail = fail+1;
    else
      pass = pass+1;
    end
  end
  
  %-----------------------------------------------------------  
  % Higher orders work over larger q domains
  % Now m = 3, 4, 5, ... 11
  qs = [-1.5, -1, -.5, -.1, .1, .5, 1, 1.5];  
  for i=1:length(qs)
    q = qs(i);
    for m=3:11
      fprintf('----------- m = %d, q = %f  -----------\n', m, q)  
      dlmfse = se_q_expansion(m,q,v);
      myse = mathieu_se(m,q,v);
      %plot(v, dlmfse,'b-')
      %hold on
      %plot(v, myse,'r.')
      %legend('dlmf','me')
      diffstd = std(dlmfse - myse);
      if (diffstd > tol)
        fprintf('Error!  m = %d, q = %5.3f, diffstd = %e\n', m, q, diffstd)
        fail = fail+1;
      else
	pass = pass+1;
      end
      %pause()
      %close all;
    end
    fprintf('--------------------------------------\n')      
  end
   
   %-----------------------------------------------------------  
  % Now even higher orders
  qs = [-10, -3, -1.5, -1, -.5, -.1, .1, .5, 1, 1.5, 3, 10];  
  for i=1:length(qs)
    q = qs(i);
    for m=12:MM
      fprintf('----------- m = %d, q = %f  -----------\n', m, q)  
      dlmfse = se_q_expansion(m,q,v);
      myse = mathieu_se(m,q,v);
      %plot(v, dlmfse,'b-')
      %hold on
      %plot(v, myse,'r.')
      %legend('dlmf','me')
      diffstd = std(dlmfse - myse);
      if (diffstd > tol)
        fprintf('Error!  m = %d, q = %5.3f, diffstd = %e\n', m, q, diffstd)
        fail = fail+1;
      else
	pass = pass+1;
      end
      %pause()
      %close all;
    end
    fprintf('--------------------------------------\n')      
  end
end

  fprintf('======================================\n')
  fprintf('At end, pass = %d, fail = %d\n', pass, fail)
  
end
