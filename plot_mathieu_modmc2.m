function plot_mathieu_modmc2()
  % This plots the modified Mathieu Mc fcns of second kind.
    
  q = 1;

  v = linspace(0,5,5000)';

  leg = {};
  
  % First do even ce fcns
  figure(1)
  for m = 0:2:6
    %fprintf('-----------------------\n')    
    y = mathieu_modmc2(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
  end  
  title('Modified Mathieu of second kind Mc2n')
  legend(leg,'Location','SouthEast')
  ylim([-6,6])
  
  leg = {};
  % Next do odd ce fcns
  figure(2)
  for m = 1:2:7
    %fprintf('-----------------------\n')
    y = mathieu_modmc2(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
  end  
  title('Modified Mathieu of second kind Mc2n+1')
  legend(leg, 'Location','SouthEast')
  ylim([-6,6])

  %==============================================
  % Compare against bessely for large u per
  % DLMF 28.20.11

  m=15;
  N = 10000;
  u = linspace(.1,5,N)';
  q = 1;
  yy = bessely(m,2*sqrt(q)*cosh(u));
  ym = mathieu_modmc2(m,q,u);
  % I need to change sign to match the Bessel fcn.
  %if (sign(ym(end)) ~= sign(yy(end)))
  %  ym = -ym;
  %end
  
  figure(4)
  plot(u,yy)
  hold on
  plot(u,ym)
  title('asymptotic behavior: modmc2 compared to bessely')
  legend('Bessely','Mathieu modmc2')
  ylim([-2,2])
  
  figure(5)
  diff = ym-yy;
  semilogy(u,abs(diff));
  title('difference')

  
end


