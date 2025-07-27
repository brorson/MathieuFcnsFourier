function plot_mathieu_modms1()
  % This plots the modified Mathieu Ms fcns of first kind.
    
  q = 1;

  v = linspace(0,5,5000);

  leg = {};
  
  % First do even Ms fcns
  figure(1)
  for m = 2:2:8
    %fprintf('-----------------------\n')    
    y = mathieu_modms1(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
  end  
  title('Modified Mathieu of first kind Ms2n')
  legend(leg)
  
  leg = {};
  % Next do odd Ms fcns
  figure(2)
  for m = 1:2:7
    %fprintf('-----------------------\n')
    y = mathieu_modms1(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
  end  
  title('Modified Mathieu of first kind Ms2n+1')
  legend(leg)

  %==============================================
  % Make plots from Guitarrez paper

  figure(3)
  m=1;
  u = linspace(0,2.5,100);
  for q = 1:3
    % I change sign to match the Guitarrez paper.
    % That said, I don't know what normalization they used,
    % so my plot doesn't exactly match theirs.
    y = -mathieu_modms1(m,q,u);
    plot(u,y)
    hold on
  end
  title('modms1 for varying q per Guitarrez')
  
  %==============================================
  % Compare against besselj for large u per
  % DLMF 28.20.11
  figure(4)
  m=5;
  N = 10000;
  u = linspace(0,5,N);
  q = 1;
  yj = besselj(m,2*sqrt(q)*cosh(u));
  ym = mathieu_modms1(m,q,u);
  % I need to change sign to match the Bessel fcn.
  if (sign(ym(end)) ~= sign(yj(end)))
    ym = -ym;
  end
  
  plot(u,yj)
  hold on
  plot(u,ym)
  title('asymptotic behavior: modms1 compared to besselj')
  legend('Besselj','Mathieu modms1')
  
  figure(5)
  diff = ym-yj;
  semilogy(u,abs(diff));
  title('difference')
  

  
end


