function plot_normalization_err()
  % This makes a contour plot of the normalization error
  % for both ce and se.

  % Parameters to vary
  ms = 1:35;
  qs = logspace(-4,4,20);
  
  % Domain
  N = 1000;
  v = linspace(0, 10, N)';

  %-----------------------------------------------------
  % 
  fprintf('Computing dot product of ce * ce\n')
  % Matrix of error values.
  errs = zeros(length(ms),length(qs));
  X = zeros(size(errs));
  Y = zeros(size(errs));  
  
  % Loop over q and m
  for i=1:length(ms)
    m = ms(i);
    fprintf('-----------  m = %d  -----------\n', m)
    for j = 1:length(qs)
      q = qs(j);

      %[y1,y1d] = mathieu_ce(m,q,v);
      %[y2,y2d] = mathieu_se(m,q,v);

      % Compute dot product
      ce1 = @(m,q,v) mathieu_ce(m,q,v)';  
      ce2 = @(m,q,v) mathieu_ce(m,q,v)';        
      f = @(v) ce1(m,q,v).*ce2(m,q,v);
      s = integral(f, -pi, pi);

      % Relative norm diff
      errs(i,j) = log10(abs(s - pi));
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  figure(2)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of normalization error -- ce and ce')


  %-----------------------------------------------------
  % 
  fprintf('Computing dot product of se * se\n')
  % Matrix of error values.
  errs = zeros(length(ms),length(qs));
  X = zeros(size(errs));
  Y = zeros(size(errs));  
  
  % Loop over q and m
  for i=1:length(ms)
    m = ms(i);
    fprintf('-----------  m = %d  -----------\n', m)
    for j = 1:length(qs)
      q = qs(j);

      %[y1,y1d] = mathieu_ce(m,q,v);
      %[y2,y2d] = mathieu_se(m,q,v);

      % Compute dot product
      se1 = @(m,q,v) mathieu_se(m,q,v)';  
      se2 = @(m,q,v) mathieu_se(m,q,v)';        
      f = @(v) se1(m,q,v).*se2(m,q,v);
      s = integral(f, -pi, pi);

      % Relative norm diff
      errs(i,j) = log10(abs(s - pi));
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  figure(3)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of normalization error -- se and se')


end
