function y = besseljd(n,x)
  % This returns the x derivative of besselj.
    
  if (floor(n) ~= n)
    error('besseljd -- order n must be integer!\n')
  end
    
  if (n == 0)
    y = -besselj(1,x);
  else
    y = (besselj(n-1,x)-besselj(n+1,x))/2;
  end

  if (n<0)
    s = (-1)^n;
    y = s*y;
  end

  return

end


