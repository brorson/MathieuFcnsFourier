function plot_mathieu_ce()
  % This plots the Mathieu ce fcns.
    
  q = 10;

  v = linspace(0,pi/2,100);

  % First do even ce fcns
  leg = {};
  figure(1)
  for m = 0:2:6
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
  for m = 1:2:7
    %fprintf('-----------------------\n')
    y = mathieu_ce(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
 end  
  title('Mathieu ce2n+1')
  legend(leg,'Location','SouthWest')
  
end


