function plot_round_trip_err()
  % This makes a contour plot of the round trip error
  % for all the Mathieu fcns.



if 0
  %========================================================
  %-----------------------------------------------------
  % ce
  MM = 100;  % Max order to test.
  v = linspace(-pi,pi*(NN-1)/NN,NN)';
  h = 5e-4;

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

      % Use sixth order deriv.  Coeffs are:
      % -1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60
      [~,fm3] = mathieu_ce(m,q,v-3*h);
      [~,fm2] = mathieu_ce(m,q,v-2*h);
      [~,fm1] = mathieu_ce(m,q,v-h);
      [~,fp1] = mathieu_ce(m,q,v+h);
      [~,fp2] = mathieu_ce(m,q,v+2*h);
      [~,fp3] = mathieu_ce(m,q,v+3*h);
      ydd = -fm3/60 + 3*fm2/20 - 3*fm1/4 + 3*fp1/4 - 3*fp2/20 + fp3/60;

      [y,yd] = mathieu_ce(m,q,v);
      a = mathieu_a(m,q);

      r = ydd/(h) + (a - 2*q*cos(2*v)).*y;
     
      stddev = std(r);
      l2norm = norm(y);
     
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
  caxis([-25 10])
  colormap default
end

if 0
  %========================================================
  %-----------------------------------------------------
  % se
  MM = 100;  % Max order to test.
  v = linspace(-pi,pi*(NN-1)/NN,NN)';
  h = 5e-4;
  
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

      % Use sixth order deriv.  Coeffs are:
      % -1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60
      [~,fm3] = mathieu_se(m,q,v-3*h);
      [~,fm2] = mathieu_se(m,q,v-2*h);
      [~,fm1] = mathieu_se(m,q,v-h);
      [~,fp1] = mathieu_se(m,q,v+h);
      [~,fp2] = mathieu_se(m,q,v+2*h);
      [~,fp3] = mathieu_se(m,q,v+3*h);
      ydd = -fm3/60 + 3*fm2/20 - 3*fm1/4 + 3*fp1/4 - 3*fp2/20 + fp3/60;

      [y,yd] = mathieu_se(m,q,v);
      a = mathieu_b(m,q);

      r = ydd/(h) + (a - 2*q*cos(2*v)).*y;
     
      stddev = std(r);
      l2norm = norm(y);
     
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
  caxis([-25 10])
  colormap default
end

if 1
  %========================================================
  %-----------------------------------------------------
  % modmc1
  MM = 50;  % Max order to test.
  NN = 100;  % Number of pts to sample in domain
  v = linspace(2,5,NN)';
  h = 5e-6;
  
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

      % Use sixth order deriv.  Coeffs are:
      % -1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60
      [~,fm3] = mathieu_modmc1(m,q,v-3*h);
      [~,fm2] = mathieu_modmc1(m,q,v-2*h);
      [~,fm1] = mathieu_modmc1(m,q,v-h);
      [~,fp1] = mathieu_modmc1(m,q,v+h);
      [~,fp2] = mathieu_modmc1(m,q,v+2*h);
      [~,fp3] = mathieu_modmc1(m,q,v+3*h);
      ydd = -fm3/60 + 3*fm2/20 - 3*fm1/4 + 3*fp1/4 - 3*fp2/20 + fp3/60;

      % Use eighth order deriv.  Coeffs are:
      % 1/280   −4/105  1/5  −4/5  0   4/5  −1/5  4/105 −1/280
      %[~,fm4] = mathieu_modmc1(m,q,v-4*h);      
      %[~,fm3] = mathieu_modmc1(m,q,v-3*h);
      %[~,fm2] = mathieu_modmc1(m,q,v-2*h);
      %[~,fm1] = mathieu_modmc1(m,q,v-h);
      %[~,fp1] = mathieu_modmc1(m,q,v+h);
      %[~,fp2] = mathieu_modmc1(m,q,v+2*h);
      %[~,fp3] = mathieu_modmc1(m,q,v+3*h);
      %[~,fp4] = mathieu_modmc1(m,q,v+4*h);      
      %ydd = fm4/280 - 4*fm3/105 + fm2/5 - 4*fm1/5 + 4*fp1/5 - fp2/5 + 4*fp3/105 - fp4/280;
  
      [y,~] = mathieu_modmc1(m,q,v);
      a = mathieu_a(m,q);

      r = ydd/(h) - (a - 2*q*cosh(2*v)).*y;
     
      stddev = std(r);
      l2norm = norm(y);
     
      % Err results
      errs(i,j) = log10(stddev/l2norm);
      X(i,j) = m;
      Y(i,j) = log10(q);
     
    end
  end

  figure(3)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of round-trip rel error -- modmc1')
  caxis([-25 10])
  colormap default
end


if 1
  %========================================================
  %-----------------------------------------------------
  % modmc2
  % Fcn is singular at v = 0;  Therefore, start comparison
  % at v = 1.
  MM = 50;   % Max order to test.
  NN = 100;  % Number of pts to sample in domain
  v = linspace(2,5,NN)';  
  h = 5e-6;
  
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

      % Use fd formula and small, fixed h to compute deriv.
      % Use sixth order deriv.  Coeffs are:
      % -1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60
      [~,fm3] = mathieu_modmc2(m,q,v-3*h);
      [~,fm2] = mathieu_modmc2(m,q,v-2*h);
      [~,fm1] = mathieu_modmc2(m,q,v-h);
      [~,fp1] = mathieu_modmc2(m,q,v+h);
      [~,fp2] = mathieu_modmc2(m,q,v+2*h);
      [~,fp3] = mathieu_modmc2(m,q,v+3*h);
      ydd = -fm3/60 + 3*fm2/20 - 3*fm1/4 + 3*fp1/4 - 3*fp2/20 + fp3/60;

      [y,~] = mathieu_modmc2(m,q,v);
      a = mathieu_a(m,q);

      r = ydd/(h) - (a - 2*q*cosh(2*v)).*y;
     
      stddev = std(r);
      l2norm = norm(y);
      
      % Err results
      errs(i,j) = log10(stddev/l2norm);
      X(i,j) = m;
      Y(i,j) = log10(q);
     
    end
  end

  figure(4)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of round-trip rel error -- modmc2')
  caxis([-25 10])
  colormap default
end


if 1
  %========================================================
  %-----------------------------------------------------
  % modms1
  MM = 50;   % Max order to test.
  NN = 100;  % Number of pts to sample in domain
  v = linspace(2,5,NN)';  
  h = 5e-6;
  
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
      
      % Use fd formula and small, fixed h to compute deriv.
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
      
      % Err results
      errs(i,j) = log10(stddev/l2norm);
      X(i,j) = m;
      Y(i,j) = log10(q);
     
    end
  end

  figure(5)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of round-trip error -- modms1')
  caxis([-25 10])
  colormap default
end


if 1
  %========================================================
  %-----------------------------------------------------
  % modms2
  MM = 50;   % Max order to test.
  NN = 100;  % Number of pts to sample in domain
  v = linspace(2,5,NN)';  
  h = 5e-6;
  
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

      % Use fd formula and small, fixed h to compute deriv.
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
     
      % Err results
      errs(i,j) = log10(stddev/l2norm);
      X(i,j) = m;
      Y(i,j) = log10(q);
     
    end
  end

  figure(6)
  contourf(X,Y,errs,-25:5:35,'ShowText','on')
  xlabel('Order m')
  ylabel('log10(q)')
  title('Log10 of round-trip error -- modms2')
  caxis([-25 10])
  colormap default
end

end
