function plot_mathieu_modms2()
  % This plots the modified Mathieu Ms fcns of second kind.
    
  q = 1;
  vmax = 4;
  v = linspace(0,vmax,1000);

  leg = {};
  
  % First do even Ms fcns
  figure(1)
  for m = 2:2:8
    %fprintf('-----------------------\n')    
    y = mathieu_modms2(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
  end  
  title('Modified Mathieu of second kind Ms2n')
  legend(leg)
  ylim([-2,2])
  
  leg = {};
  % Next do odd Ms fcns
  figure(2)
  for m = 1:2:7
    %fprintf('-----------------------\n')
    y = mathieu_modms2(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
  end  
  title('Modified Mathieu of second kind Ms2n+1')
  legend(leg)
  ylim([-2,2])  

  %==============================================
  % Compare against besselj for large u per
  % DLMF 28.20.11
  figure(4)
  m=7;
  N = 1000;
  u = linspace(0,vmax,N);
  q = 1;
  yy = bessely(m,2*sqrt(q)*cosh(u));
  ym = mathieu_modms2(m,q,u);
  % I need to change sign to match the Bessel fcn.
  %if (sign(ym(end)) ~= sign(yy(end)))
  %  ym = -ym;
  %end
  
  plot(u,yy)
  hold on
  plot(u,ym)
  title('asymptotic behavior: modms2 compared to bessely')
  legend('Bessely','Mathieu modms2')
  ylim([-6,6])
  
  figure(5)
  diff = ym-yy;
  semilogy(u,abs(diff));
  title('difference')

  
end


