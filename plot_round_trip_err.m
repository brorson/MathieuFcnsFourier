function plot_round_trip_err()
  % This makes a contour plot of the round trip error
  % for fcns ce and se.

  % Domain
  % 3000 seems to be the sweet spot using 8th order deriv.
  NN = 3000;  % More pts => smaller h, more accurate test.
              % But at some point the accuracy turns around ...
	      % not sure why.
  MM = 30;  % Max order to test.

if 1
  %========================================================
  %-----------------------------------------------------
  % ce
  v = linspace(-pi,pi*(NN-1)/NN,NN)';
  h = v(2)-v(1);

  % Parameters to vary
  ms = 0:MM;  % ce orders start at 0
  qs = logspace(-4,4,30);
  
  fprintf('Computing round trip error for ce\n')
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

      [y,yd] = mathieu_ce(m,q,v);
      a = mathieu_a(m,q);

     ydd = fd_deriv(yd);
     r = ydd/(h) + (a - 2*q*cos(2*v)).*y;
     
     % I have some sort of bug in taking the deriv, so I need
     % to truncate the residual vector.  Must fix this later.
     stddev = std(r(5:(end-4)));
     l2norm = norm(y(5:end-4));
     
     % Err results
     errs(i,j) = log10(stddev/l2norm);
     X(i,j) = m;
     Y(i,j) = log10(q);
     
    end
  end

  figure(1)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of round-trip error -- ce')
  caxis([-20 5])
  colormap default
end

if 1
  %========================================================
  %-----------------------------------------------------
  % se
  v = linspace(-pi,pi*(NN-1)/NN,NN)';
  h = v(2)-v(1);

  % Parameters to vary
  ms = 1:MM;  % se orders start at 1
  qs = logspace(-4,4,30);
  
  fprintf('Computing round trip error for se\n')
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

      [y,yd] = mathieu_se(m,q,v);
      a = mathieu_b(m,q);

     ydd = fd_deriv(yd);
     r = ydd/(h) + (a - 2*q*cos(2*v)).*y;
     
     % I have some sort of bug in taking the deriv, so I need
     % to truncate the residual vector.  Must fix this later.
     stddev = std(r(5:(end-4)));
     l2norm = norm(y(5:end-4));
     
     % Err results
     errs(i,j) = log10(stddev/l2norm);
     X(i,j) = m;
     Y(i,j) = log10(q);
     
    end
  end

  figure(2)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of round-trip error -- se')
  caxis([-20 5])
  colormap default
end


  %========================================================
  %-----------------------------------------------------
  % modmc1
  v = linspace(0,5,NN)';
  h = v(2)-v(1);

  % Parameters to vary
  ms = 1:MM;  % se orders start at 1
  qs = logspace(-4,4,10);
  
  fprintf('Computing round trip error for modmc1\n')
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

      [y,yd] = mathieu_modmc1(m,q,v);
      a = mathieu_a(m,q);

     ydd = fd_deriv(yd);
     r = ydd/(h) - (a - 2*q*cosh(2*v)).*y;
     
     % I have some sort of bug in taking the deriv, so I need
     % to truncate the residual vector.  Must fix this later.
     stddev = std(r(5:(end-4)));
     
     % Err results
     errs(i,j) = log10(stddev);
     X(i,j) = m;
     Y(i,j) = log10(q);
     
    end
  end

  figure(3)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of round-trip error -- modmc1')
  caxis([-20 5])
  colormap default

  
  %========================================================
  %-----------------------------------------------------
  % modmc2
  v = linspace(0,5,NN)';
  h = v(2)-v(1);

  % Parameters to vary
  ms = 1:MM;  % se orders start at 1
  qs = logspace(-4,4,10);
  
  fprintf('Computing round trip error for modmc2\n')
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

      [y,yd] = mathieu_modmc2(m,q,v);
      a = mathieu_a(m,q);

     ydd = fd_deriv(yd);
     r = ydd/(h) - (a - 2*q*cosh(2*v)).*y;
     
     % I have some sort of bug in taking the deriv, so I need
     % to truncate the residual vector.  Must fix this later.
     stddev = std(r(5:(end-4)));
     
     % Err results
     errs(i,j) = log10(stddev);
     X(i,j) = m;
     Y(i,j) = log10(q);
     
    end
  end

  figure(4)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of round-trip error -- modmc2')
  caxis([-20 5])
  colormap default


  
  %========================================================
  %-----------------------------------------------------
  % modms1
  v = linspace(0,5,NN)';
  h = v(2)-v(1);

  % Parameters to vary
  ms = 1:MM;  % se orders start at 1
  qs = logspace(-4,4,10);
  
  fprintf('Computing round trip error for modms1\n')
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

      [y,yd] = mathieu_modms1(m,q,v);
      a = mathieu_b(m,q);

     ydd = fd_deriv(yd);
     r = ydd/(h) - (a - 2*q*cosh(2*v)).*y;
     
     % I have some sort of bug in taking the deriv, so I need
     % to truncate the residual vector.  Must fix this later.
     stddev = std(r(5:(end-4)));
     
     % Err results
     errs(i,j) = log10(stddev);
     X(i,j) = m;
     Y(i,j) = log10(q);
     
    end
  end

  figure(5)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of round-trip error -- modms1')
  caxis([-20 5])
  colormap default


  %========================================================
  %-----------------------------------------------------
  % modms2
  v = linspace(0,5,NN)';
  h = v(2)-v(1);

  % Parameters to vary
  ms = 1:MM;  % se orders start at 1
  qs = logspace(-4,4,10);
  
  fprintf('Computing round trip error for modms2\n')
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

      [y,yd] = mathieu_modms2(m,q,v);
      a = mathieu_b(m,q);

     ydd = fd_deriv(yd);
     r = ydd/(h) - (a - 2*q*cosh(2*v)).*y;
     
     % I have some sort of bug in taking the deriv, so I need
     % to truncate the residual vector.  Must fix this later.
     stddev = std(r(5:(end-4)));
     
     % Err results
     errs(i,j) = log10(stddev);
     X(i,j) = m;
     Y(i,j) = log10(q);
     
    end
  end

  figure(6)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of round-trip error -- modms2')
  caxis([-20 5])
  colormap default


end
