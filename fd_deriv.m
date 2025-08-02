function y = fd_deriv(x)
  % This computes the finite difference deriv of input vector
  % x.  It returns the deriv in vector y which has the same
  % length as x.  It assumes cyclic boundary conds.
  % It also assumes the right end point is to the left of the
  % first point.  Therefore, linspace(-pi,pi,N) must have the
  % right point trimmed before calling fd_deriv.
    
  ord = 8;

  if (ord == 4)
    % Fourth order expression.  Coeffs are 
    % 1/12, -2/3, 0, 2/3, -1/12
    y = zeros(size(x));
    y(1) = x(end-1)/12 - 2*x(end)/3 + 2*x(2)/3 - x(3)/12;
    y(2) = x(end)/12 - 2*x(1)/3 + 2*x(3)/3 - x(4)/12;      
    
    y(3:end-2) = x(1:end-4)/12 - 2*x(2:end-3)/3 + 2*x(4:end-1)/3 - x(5:end)/12;
    
    y(end-1) = x(end-3)/12 - 2*x(end-2)/3 + 2*x(end)/3 - x(1)/12;      
    y(end) = x(end-2)/12 - 2*x(end-1)/3 + 2*x(1)/3 - x(2)/12;            
    
  elseif (ord == 6)
    % Sixth order.  Coeffs are:
    % -1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60
    y = zeros(size(x));
    y(1) = -x(end-2)/60 + 3*x(end-1)/20 - 3*x(end)/4 + 3*x(2)/4 - 3*x(3)/20 + x(4)/60;
    y(2) = -x(end-1)/60 + 3*x(end)/20   - 3*x(1)/4   + 3*x(3)/4 - 3*x(4)/20 + x(5)/60;
    y(3) = -x(end)/60   + 3*x(1)/20     - 3*x(2)/4   + 3*x(4)/4 - 3*x(5)/20 + x(6)/60;    
    
    y(4:end-3) = -x(1:end-6)/60 + 3*x(2:end-5)/20 - 3*x(3:end-4)/4 + 3*x(5:end-2)/4 - 3*x(6:end-1)/20 + x(7:end)/60;        
    
    y(end-2) = -x(end-5)/60 + 3*x(end-4)/20 - 3*x(end-3)/4 + 3*x(end-1)/4 - 3*x(end)/20 + x(1)/60;        
    y(end-1) = -x(end-4)/60 + 3*x(end-3)/20 - 3*x(end-2)/4 + 3*x(end)/4 - 3*x(1)/20 + x(2)/60;        
    y(end) = -x(end-3)/60 + 3*x(end-2)/20 - 3*x(end-1)/4 + 3*x(1)/4 - 3*x(2)/20 + x(3)/60;        

  elseif(ord == 8)
    % 8th order.  Coeffs = 
    % 1/280 	−4/105 	1/5 	−4/5 	0 	4/5 	−1/5 	4/105 	−1/280
    y = zeros(size(x));
    y(1) = x(end-3)/280 - 4*x(end-2)/105 + x(end-1)/5 - 4*x(end)/5 ...
	   + 4*x(2)/5 - x(3)/5 + 4*x(4)/105 - x(5)/280;
    y(2) = x(end-2)/280 - 4*x(end-1)/105 + x(end)/5 - 4*x(1)/5 ...
	   + 4*x(3)/5 - x(4)/5 + 4*x(5)/105 - x(6)/280;
    y(3) = x(end-1)/280 - 4*x(end)/105 + x(1)/5 - 4*x(2)/5 ...
	   + 4*x(4)/5 - x(5)/5 + 4*x(6)/105 - x(7)/280;
    y(4) = x(end)/280 - 4*x(1)/105 + x(2)/5 - 4*x(3)/5 ...
	   + 4*x(5)/5 - x(6)/5 + 4*x(7)/105 - x(8)/280;

    y(5:end-4) = x(1:end-8)/280 - 4*x(2:end-7)/105 + x(3:end-6)/5 - 4*x(4:end-5)/5 ...
	   + 4*x(6:end-3)/5 - x(7:end-2)/5 + 4*x(8:end-1)/105 - x(9:end)/280;

    y(end-3) = x(end-7)/280 - 4*x(end-6)/105 + x(end-5)/5 - 4*x(end-4)/5 ...
	   + 4*x(end-2)/5 - x(end-1)/5 + 4*x(end)/105 - x(1)/280;
    y(end-2) = x(end-6)/280 - 4*x(end-5)/105 + x(end-4)/5 - 4*x(end-3)/5 ...
	   + 4*x(end-1)/5 - x(end)/5 + 4*x(1)/105 - x(2)/280;
    y(end-1) = x(end-5)/280 - 4*x(end-4)/105 + x(end-3)/5 - 4*x(end-2)/5 ...
	   + 4*x(end)/5 - x(1)/5 + 4*x(2)/105 - x(3)/280;
    y(end) = x(end-4)/280 - 4*x(end-3)/105 + x(end-2)/5 - 4*x(end-1)/5 ...
	   + 4*x(1)/5 - x(2)/5 + 4*x(3)/105 - x(4)/280;
	
  else
    error('Bad order chosen for fd_deriv!')
  end
  
end


      