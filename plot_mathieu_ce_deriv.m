function plot_mathieu_ce_deriv()
  % This plots the Mathieu ce_deriv fcns.
    
  q = 1;
  
  fprintf('Plotting Fourier series ce deriv for q = %f\n', q)
  
  v = linspace(-pi,pi,100)';

  % First do even ce deriv fcns
  leg = {};
  figure(1)
  for m = 0:2:6
    %fprintf('-----------------------\n')    
    [~,y] = mathieu_ce(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
   end  
  title('Mathieu ce2n deriv')
  legend(leg,'Location','SouthWest')
  
  % Next do odd ce deriv fcns
  leg = {};
  figure(2)
  for m = 1:2:7
    %fprintf('-----------------------\n')
    [~,y] = mathieu_ce(m,q,v);
    plot(v,y)
    hold on
    ss = ['m = ',num2str(m)];
    leg = [leg, ss];
 end  
  title('Mathieu ce2n+1 deriv')
  legend(leg,'Location','SouthWest')
  
end


