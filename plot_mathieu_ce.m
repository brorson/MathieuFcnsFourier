function plot_mathieu_ce()
  % This plots the Mathieu ce fcns.
    
  q = 1;
  
  fprintf('Plotting Fourier series ce for q = %f\n', q)
  
  v = linspace(-pi,pi,100);

  % First do even ce fcns
  leg = {};
  figure(1)
  for m = 0:2:4
    %fprintf('-----------------------\n')    
    y = mathieu_ce(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
   end  
  title('Mathieu ce2n')
  legend(leg,'Location','SouthWest')
  
  % Next do odd ce fcns
  leg = {};
  figure(2)
  for m = 1:2:5
    %fprintf('-----------------------\n')
    y = mathieu_ce(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
 end  
  title('Mathieu ce2n+1')
  legend(leg,'Location','SouthWest')

  %=========================================================
    %==============================================
  % Compare against high-order FD approx to original eq.
  % DLMF 28.20.11
  m=20;
  N = 1000;
  u = linspace(0,5,N)';
  h = u(2)-u(1);
  q = 100;

  [y,yd] = mathieu_ce(m,q,u);
  a = mathieu_a(m,q);

  ydd = fd_deriv(yd);
  r = ydd/(h) + (a - 2*q*cos(2*u)).*y;
     
  plt_range = 5:N-4;
  
  figure(6)
  plot(u,y)
  hold on
  plot(u,yd)
  title('ce and ced')
  
  figure(7)
  plot(u(plt_range), r(plt_range))
  title('Deviation from finite diff approx')

  err = std(r(plt_range));
  l2norm = norm(r(plt_range));
  fprintf('Order m = %d, parameter q = %f, abs err = %e, norm err = %e\n', ...
	  m, q, err, err/l2norm)
  
end


