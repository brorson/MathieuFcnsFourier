function plot_orthogonality_err()
  % This makes a contour plot of the orthogonality error
  % between ce and se.

  % Parameters to vary
  ms = 1:50;
  qs = logspace(-5,5,50);
  
  c = 5;
  
  % Domain
  N = 1000;

  
  %-----------------------------------------------------
  % 
  fprintf('Computing dot product of ce m* se m+%d\n',c)
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
      ce = @(p,q,v) mathieu_ce(p,q,v)';  
      se = @(p,q,v) mathieu_se(p+c,q,v)';        
      f = @(v) ce(m,q,v).*se(m,q,v);
      s = integral(f, -pi, pi);

      % Relative norm diff
      errs(i,j) = log10(abs(s));
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  figure(1)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Base order m')
  ylabel('log10(q)')
  title(['Log10 of orthogonality error -- ce m and se m+',num2str(c)])

  %-----------------------------------------------------
  % 
  fprintf('Computing dot product of ce m * ce m + %d\n', c)
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

      % Compute dot product
      ce1 = @(p,q,v) mathieu_ce(p,q,v)';  
      ce2 = @(p,q,v) mathieu_ce(p+c,q,v)';        
      f = @(v) ce1(m,q,v).*ce2(m,q,v);
      s = integral(f, -pi, pi);

      % Relative norm diff
      errs(i,j) = log10(abs(s));
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  % Make a nice-looking plot by changing -inf to
  % -25
  idx = find(errs<-25);
  errs(idx) = -24.999;
  
  figure(2)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Base order m')
  ylabel('log10(q)')
  title(['Log10 of orthogonality error -- ce m and ce m+',num2str(c)])


  %-----------------------------------------------------
  % 
  fprintf('Computing dot product of se m * se m + %d\n', c)
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

      % Compute dot product
      se1 = @(p,q,v) mathieu_se(p,q,v)';  
      se2 = @(p,q,v) mathieu_se(p+c,q,v)';        
      f = @(v) se1(m,q,v).*se2(m,q,v);
      s = integral(f, -pi, pi);

      % Relative norm diff
      errs(i,j) = log10(abs(s));
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  % Make a nice-looking plot by changing -inf to
  % -25
  idx = find(errs<-25);
  errs(idx) = -24.999;
  
  figure(3)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Base order m')
  ylabel('log10(q)')
  title(['Log10 of orthogonality error -- se m and se m+',num2str(c)])

end
