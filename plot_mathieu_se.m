function plot_mathieu_se()
  % This plots the Mathieu se fcns.
    
  q = 1;

  v = linspace(0,pi/2,100);

  % First do even se fcns
  leg = {};
  figure(1)
  for m = 2:2:8
    %fprintf('-----------------------\n')    
    y = mathieu_se(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
  end  
  title('Mathieu se2n')
  legend(leg,'Location','SouthWest')

  % Next do odd se fcns
  leg = {};
  figure(2)
  for m = 1:2:7
    %fprintf('-----------------------\n')
    y = mathieu_se(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
  end  
  title('Mathieu se2n+1')
  legend(leg,'Location','SouthWest')
  
end


