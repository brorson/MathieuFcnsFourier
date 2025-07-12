function playpen()

if 0
  % This is a test that the eigenvalue is a constant.
  close all;

  fprintf('Testing ce to verify eigenvalue is const\n')
  m = 1;
  q = .1;
  u = linspace(-2*pi+1,2*pi,5000);
  du = u(2)-u(1);

  D2S = (mathieu_ce(m,q,u(3:end))-2*mathieu_ce(m,q,u(2:end-1))+mathieu_ce(m,q,u(1:end-2)))/(du*du);
  S = mathieu_ce(m,q,u(2:end-1));

  a = -D2S./S + 2*q*cos(2*u(2:end-1));
  fprintf('Minmax = %e\n', max(a)-min(a))

  figure(1)
  semilogy(u(2:end-1), abs(a))
  title('ce -- Eigenvalue a (semilogy)')

  figure(2)
  plot(u(2:end-1), a)
  title('ce -- Eigenvalue a (plot)')

end


if 0
  % This is a test that the eigenvalue is a constant.
  fprintf('Testing se to verify eigenvalue is const\n')
  
  m = 2;
  q = .1;
  u = linspace(-2*pi+1,2*pi,5000);
  du = u(2)-u(1);

  D2S = (mathieu_se(m,q,u(3:end))-2*mathieu_se(m,q,u(2:end-1))+mathieu_se(m,q,u(1:end-2)))/(du*du);
  S = mathieu_se(m,q,u(2:end-1));

  b = -D2S./S + 2*q*cos(2*u(2:end-1));
  fprintf('Minmax = %e\n', max(b)-min(b))

  figure(11)
  semilogy(u(2:end-1), abs(b))
  title('se -- Eigenvalue b (semilogy)')

  figure(12)
  plot(u(2:end-1), b)
  title('se -- Eigenvalue b (plot)')
  %ylim([-10,10])

end


if 1
  % This is a test that the eigenvalue is a constant.
  fprintf('Testing modce1 to verify eigenvalue is const\n')

  m = 25;
  q = 1;
  u = linspace(0,10,50000);
  du = u(2)-u(1);

  D2R = (mathieu_modce1(m,q,u(3:end))-2*mathieu_modce1(m,q,u(2:end-1))+mathieu_modce1(m,q,u(1:end-2)))/(du*du);
  R = mathieu_modce1(m,q,u(2:end-1));

  figure(100)
  plot(u(2:end-1), D2R)
  title('2nd deriv of modce1')
  
  figure(101)
  plot(u(2:end-1), R)
  title('modce1')
  
  a = D2R./R + 2*q*cosh(2*u(2:end-1));
  fprintf('Minmax = %e\n', max(a)-min(a))

  figure(3)
  semilogy(u(2:end-1), abs(a))
  title('modce1 -- Eigenvalue a (semilogy)')

  figure(4)
  plot(u(2:end-1), a)
  title('modce1 -- Eigenvalue a (plot)')

end


if 0
  %======================================================
  % Test diff eq for mathieu_ce

  % close all;
  fprintf('Test diff eq for mathieu_ce\n')
  
  m = 0;
  q = 1;
  v = linspace(-pi+1,pi,500);
  dv = v(2)-v(1);

if 0
  D2S = (mathieu_ce(m,q,v(3:end))-2*mathieu_ce(m,q,v(2:end-1))+mathieu_ce(m,q,v(1:end-2)))/(dv*dv);
  S = mathieu_ce(m,q,v(2:end-1));

  a = -D2S./S + 2*q*cos(2*v(2:end-1));
  fprintf('Minmax = %e\n', max(a)-min(a))

  figure(1)
  semilogy(v(2:end-1), abs(a))
  title('Eigenvalue a (semilogy)')

  figure(2)
  plot(v(2:end-1), a)
  title('Eigenvalue a (plot)')

  figure(3)
  plot(v(2:end-1),D2S)
  hold on
  plot(v(2:end-1),S)


  v1 = linspace(-pi+1,pi,500);
  dv1 = v1(2)-v1(1);

  D2S1 = (mathieu_ce(m,q,v1(3:end))-2*mathieu_ce(m,q,v1(2:end-1))+mathieu_ce(m,q,v1(1:end-2)))/(dv1*dv1);
  S1 = mathieu_ce(m,q,v1(2:end-1));
  
  figure(4)
  plot(v1(2:end-1),D2S1+1.3,'bo')
  hold on
  plot(v1(2:end-1),-6*S1+3.8,'ro')
  plot(v1(2:end-1), 2*q*cos(2*v1(2:end-1)), '-');
end
  
  
  figure(5)
  se1 = mathieu_se(m,q,v(3:end));
  se2 = mathieu_se(m,q,v(1:end-2));
  plot(v(3:end),se1)
  hold on
  plot(v(1:end-2),se2)
  title('Compare eval of Mathieu at two different point sets')
  
  
  figure(6)
  se1 = mathieu_se(m,q,v(3:end));
  se2 = mathieu_se(m,q,v(1:end-2));
  D1 = (se1-se2)/(2*dv);
  plot(v(2:end-1),D1)
  hold on
  title('Mathieu first deriv')

  figure(7)
  se1 = mathieu_se(m,q,v(3:end));
  se2 = mathieu_se(m,q,v(2:end-1));
  se3 = mathieu_se(m,q,v(1:end-2));
  D2 = (se1-2*se2+se3)/(dv*dv);
  plot(v(2:end-1),D2)
  hold on
  title('Mathieu second deriv')
  
end


if 0
  %======================================================
  % Test diff eq for cos(v)

  close all;

  m = 0;
  q = 1;
  v = linspace(-2,2*pi,500);
  dv = v(2)-v(1);

  D2S = ( cos(v(3:end)) - 2*cos(v(2:end-1)) + cos(v(1:end-2)) )/(dv*dv);
  S = cos(v(2:end-1));

  a = D2S./S;
  fprintf('Minmax = %e\n', max(a)-min(a))

  figure(1)
  semilogy(v(2:end-1), abs(a))

  figure(2)
  plot(v(2:end-1), a)
  
end

end
