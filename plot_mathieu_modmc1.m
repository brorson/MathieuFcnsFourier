function plot_mathieu_modmc1()
  % This plots the modified Mathieu Mc fcns of first kind.
    
  q = 30;

  v = linspace(0,5,1000)';

  leg = {};
  
  % First do even ce fcns
  figure(1)
  for m = 0:2:6
    %fprintf('-----------------------\n')    
    y = mathieu_modmc1(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
  end  
  title('Modified Mathieu of first kind Mc2n')
  legend(leg,'Location','SouthEast')

  leg = {};
  % Next do odd ce fcns
  figure(2)
  for m = 1:2:7
    %fprintf('-----------------------\n')
    y = mathieu_modmc1(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
  end  
  title('Modified Mathieu of first kind Mc2n+1')
  legend(leg, 'Location','SouthEast')

  %==============================================
  % Make plots of Guitarrez paper

  figure(3)
  m=1;
  u = linspace(0,2.5,100)';
  for q = 1:3
    % I change sign to match the Guitarrez paper 
    y = -mathieu_modmc1(m,q,u);
    plot(u,y)
    hold on
  end
  title('modmc1 for varying q per Guitarrez')

  %==============================================
  % Compare against besselj for large u per
  % DLMF 28.20.11
  m=5;
  N = 10000;
  u = linspace(0,5,N)';
  q = 30;
  yj = besselj(m,2*sqrt(q)*cosh(u));
  ym = mathieu_modmc1(m,q,u);
  % I need to change sign to match the Bessel fcn.
  %if (sign(ym(end)) ~= sign(yj(end)))
  %  ym = -ym;
  %end

  figure(4)
  plot(u,yj)
  hold on
  plot(u,ym)
  title('asymptotic behavior: modmc1 compared to besselj')
  legend('Besselj','Mathieu modmc1')
  
  figure(5)
  diff = ym-yj;
  semilogy(u,abs(diff));
  title('difference')
    
end


