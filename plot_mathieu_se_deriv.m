function plot_mathieu_se_deriv()
  % This plots the Mathieu se_deriv fcns.
    
  q = -1;

  v = linspace(-pi,pi,100);

  % First do even se_deriv fcns
  leg = {};
  figure(1)
  for m = 2:2:2
    %fprintf('-----------------------\n')    
    y = mathieu_se_deriv(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
  end  
  title('Mathieu se2n deriv')
  legend(leg,'Location','SouthWest')

  % Next do odd se fcns
  leg = {};
  figure(2)
  for m = 1:2:1
    %fprintf('-----------------------\n')
    y = mathieu_se_deriv(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
  end  
  title('Mathieu se2n+1 deriv')
  legend(leg,'Location','SouthWest')
  
end


