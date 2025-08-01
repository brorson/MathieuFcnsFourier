function plot_wronskian_err()
  % This makes a contour plot of the error

  % Parameters to vary
  ms = 0:20;
  qs = logspace(-3,3,20);
  
  % Domain
  N = 1000;
  v = linspace(0, 10, N)';

  % True value of Wronskian  
  wtrue = 2/pi;  

  %-----------------------------------------------------
  % Mc1 & Mc2
  fprintf('Computing Wronskian of Mc1 & Mc2\n')
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

      [y1,y1d] = mathieu_modmc1(m,q,v);
      [y2,y2d] = mathieu_modmc2(m,q,v);

      % Compute Wronskian
      w = y1.*y2d - y1d.*y2;
      
      % Relative norm diff
      errs(i,j) = log10(norm(w-wtrue)/N);
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  figure(1)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of Wronskian error -- Mc1 Mc2')

  %-----------------------------------------------------
  % Ms1 & Ms2
  fprintf('Computing Wronskian of Ms1 & Ms2\n')

  % Ms supports orders 1, 2, 3, ...
  ms = 1:21;

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

      [y1,y1d] = mathieu_modms1(m,q,v);
      [y2,y2d] = mathieu_modms2(m,q,v);

      % Compute Wronskian
      w = y1.*y2d - y1d.*y2;
      
      % Relative norm diff
      errs(i,j) = log10(norm(w-wtrue)/N);
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  figure(2)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of Wronskian error -- ms1 ms2')

  %-----------------------------------------------------
  % Mc1 & Ms2
  fprintf('Computing Wronskian of Mc1 & Ms2\n')

  % Ms supports orders 1, 2, 3, ...
  ms = 1:21;

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

      [y1,y1d] = mathieu_modmc1(m,q,v);
      [y2,y2d] = mathieu_modms2(m,q,v);

      % Compute Wronskian
      w = y1.*y2d - y1d.*y2;
      
      % Relative norm diff
      errs(i,j) = log10(norm(w-wtrue)/N);
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  figure(3)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of Wronskian error -- mc1 ms2')
  %-----------------------------------------------------
  % Ms1 & Mc2
  fprintf('Computing Wronskian of Ms1 & Mc2\n')

  % Ms supports orders 1, 2, 3, ...
  ms = 1:21;

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

      [y1,y1d] = mathieu_modms1(m,q,v);
      [y2,y2d] = mathieu_modmc2(m,q,v);

      % Compute Wronskian
      w = y1.*y2d - y1d.*y2;
      
      % Relative norm diff
      errs(i,j) = log10(norm(w-wtrue)/N);
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  figure(4)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of Wronskian error -- ms1 mc2')

end
