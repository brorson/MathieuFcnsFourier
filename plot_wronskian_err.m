function plot_wronskian_err()
  % This makes a contour plot of the error

  % Parameters to vary
  ms = 1:50;
  qs = logspace(-4,4,30);
  
  % Domain
  N = 100;
  v = linspace(5, 10, N)';

  % True value of Wronskian  
  wtrue = 2/pi;  

if 1
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
      
      % Err results
      errs(i,j) = log10(std(w-wtrue));
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  figure(1)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of Wronskian error -- Mc1 Mc2')
  caxis([-20 10])
  colormap default
end


if 1
  %-----------------------------------------------------
  % Ms1 & Ms2
  fprintf('Computing Wronskian of Ms1 & Ms2\n')

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
      
      % Err results
      errs(i,j) = log10(std(w-wtrue));
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  figure(2)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of Wronskian error -- ms1 ms2')
  caxis([-20 10])
  colormap default
end


if 1
  %-----------------------------------------------------
  % Mc1 & Ms2
  fprintf('Computing Wronskian of Mc1 & Ms2\n')

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
      
      % Err results
      errs(i,j) = log10(std(w-wtrue));
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  figure(3)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of Wronskian error -- mc1 ms2')
  caxis([-20 10])
  colormap default
end


if 1
  %-----------------------------------------------------
  % Ms1 & Mc2
  fprintf('Computing Wronskian of Ms1 & Mc2\n')

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
      
      % Err results
      errs(i,j) = log10(std(w-wtrue));
      X(i,j) = m;
      Y(i,j) = log10(q);
      
    end
  end

  figure(4)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of Wronskian error -- ms1 mc2')
  caxis([-20 10])
  colormap default
end


end
