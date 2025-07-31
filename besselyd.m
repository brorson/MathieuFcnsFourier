function y = besselyd(n,x)
  % This returns the x derivative of bessely.
    
  if (floor(n) ~= n)
    error('besselyd -- order n must be integer!\n')
  end
    
  if (n == 0)
    y = -bessely(1,x);
    return
  else
    y = (bessely(n-1,x)-bessely(n+1,x))/2;
    return
  end

end
