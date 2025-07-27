function plot_mathieu_ce_deriv()
  % This plots the Mathieu ce_deriv fcns.
    
  q = 100;
  
  fprintf('Plotting Fourier series ce for q = %f\n', q)
  
  v = linspace(-pi,pi,100);

  % First do even ce fcns
  leg = {};
  figure(1)
  for m = 0:2:0
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
  for m = 1:2:1
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


