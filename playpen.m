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


if 0
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

if 0
  % THis reproduces the sign-flip bug.
  fprintf('Sign flip bug for coeffs_oo\n')
  q = linspace(-100,100,100);
  m = 3;
  N = 25;
  for i=1:100;
    A = mathieu_coeffs_oo(N,q(i),m);
    y(i) = A(3);
  end;
  close all;
  plot(q,y)
end


if 1

  fprintf('Sign flip bug for coeffs_ee\n')
  q = linspace(-100,100,100);
  y0 = zeros(size(q));
  y2 = zeros(size(q));
  y4 = zeros(size(q));
  m = 2;
  N = 25;
  for i=1:100;
    A = mathieu_coeffs_ee(N,q(i),m);
    y0(i) = A(1);
    y2(i) = A(2);
    y4(i) = A(3);
    y6(i) = A(4);
  end;
  close all;
  plot(q,y0)
  hold on
  plot(q,y2)
  plot(q,y4)
  plot(q,y6)
  
  m=0;
  for i=1:100;
    A = mathieu_coeffs_ee(N,q(i),m);
    y0(i) = A(1);
    y2(i) = A(2);
    y4(i) = A(3);
    y6(i) = A(4);
  end;
  plot(q,y0)
  plot(q,y2)
  plot(q,y4)
  plot(q,y6)
  legend('m2,1','m2,2','m2,3','m2,4','m0,1','m0,2','m0,3','m0,4')
  
end

if 0
  fprintf('Sign flip bug for coeffs_eo\n')
  q = linspace(-100,100,100);
  m = 3;
  N = 25;
  for i=1:100;
    A = mathieu_coeffs_eo(N,q(i),m);
    y(i) = A(3);
  end;
  close all;
  plot(q,y)
end


if 0
  fprintf('Sign flip bug for coeffs_oe\n')
  q = linspace(-100,100,100);
  m = 4;
  N = 25;
  for i=1:100;
    A = mathieu_coeffs_oe(N,q(i),m);
    y(i) = A(3);
  end;
  close all;
  plot(q,y)
end



% Failing case
%v = linspace(-pi,pi,1000);
%ce = mathieu_ce(1,1000,v);
%close all; plot(v,ce)
%se = mathieu_se(1,-1000,pi/2-v);
%hold on; plot(v,se)
 
% Another failing case 
v = linspace(-pi,pi,1000);
ce1 = mathieu_ce(2,-100,v);
close all; plot(v,ce1)
ce2 = mathieu_ce(2,100,pi/2-v);  % Incorrect -> is copy of ce1 for q=-1000
hold on; plot(v,ce2)
ce3 = mathieu_ce(2,100,v);
hold on; plot(v,ce3)
legend('ce1','ce2','ce3')

%-------------------------------------------------------
close all; 
figure(1)
v = linspace(-pi,pi,1000);
ce2p = mathieu_ce(2,630,v);
plot(v,ce2p)
ce2m = mathieu_ce(2,-630,v);
hold on; plot(v,ce2m)
legend('+630','-630')
title('ce m=2, q=630')
fprintf('ce(v=0)=%e\n', mathieu_ce(2,630,0))
fprintf('ce(v=0)=%e\n', mathieu_ce(2,-630,0))

figure(2)
ce2p = mathieu_ce(2,200,v);
plot(v,ce2p)
ce2m = mathieu_ce(2,-200,v);
hold on; plot(v,ce2m)
legend('200','-200')
title('ce m=2, q=200')
fprintf('ce(v=0)=%e\n', mathieu_ce(2,200,0))
fprintf('ce(v=0)=%e\n', mathieu_ce(2,-200,0))


figure(2)
ce3p = mathieu_ce(3,100,v);
plot(v,ce3p)
ce3m = mathieu_ce(3,-100,v);
hold on; plot(v,ce3m)
legend('+100','-100')
title('ce m=3, q=100')

figure(5)
ce4p = mathieu_ce(4,100,v);
plot(v,ce4p)
ce4m = mathieu_ce(4,-100,v);
hold on; plot(v,ce4m)
legend('+100','-100')
title('ce m=4, q=100')
 
figure(6)
ce5p = mathieu_ce(5,100,v);
plot(v,ce5p)
ce5m = mathieu_ce(5,-100,v);
hold on; plot(v,ce5m)
legend('+100','-100')
title('ce m=5, q=100')
 

figure(7)
se2p = mathieu_se(2,100,v);
plot(v,se2p)
se2m = mathieu_se(2,-100,v);
hold on; plot(v,se2m)
legend('+100','-100')
title('se m=2, q=100')

figure(8)
se3p = mathieu_se(3,100,v);
plot(v,se3p)
se3m = mathieu_se(3,-100,v);
hold on; plot(v,se3m)
legend('+100','-100')
title('se m=3, q=100')

figure(9)
se4p = mathieu_se(4,100,v);
plot(v,se4p)
se4m = mathieu_se(4,-100,v);
hold on; plot(v,se4m)
legend('+100','-100')
title('se m=4, q=100')

figure(10)
se5p = mathieu_se(5,100,v);
plot(v,se5p)
se5m = mathieu_se(5,-100,v);
hold on; plot(v,se5m)
legend('+100','-100')
title('se m=5, q=100')

%figure(3)
%ce2p = mathieu_ce(2,1000,v);
%plot(v,ce2p)
%ce2m = mathieu_ce(2,-1000,v);
%hold on; plot(v,ce2m)
%legend('+1000','-1000')
%title('ce m=2, q=1000')
 
%figure(4)
%ce3p = mathieu_ce(3,1000,v);
%plot(v,ce3p)
%ce3m = mathieu_ce(3,-1000,v);
%hold on; plot(v,ce3m)
%legend('+1000','-1000')
%title('ce m=3, q=1000')
 


%figure(9)
%se2p = mathieu_se(2,1000,v);
%plot(v,se2p)
%se2m = mathieu_se(2,-1000,v);
%hold on; plot(v,se2m)
%legend('+1000','-1000')
%title('se m=2, q=1000')
 
%figure(10)
%se3p = mathieu_se(3,1000,v);
%plot(v,se3p)
%se3m = mathieu_se(3,-1000,v);
%hold on; plot(v,se3m)
%legend('+1000','-1000')
%title('se m=3, q=1000')
 


%-------------------------------------------------------
%close all; 
 
% Another failing case 
%v = linspace(-pi,pi,1000);
%ce1 = mathieu_ce(2,-1,v);
%close all; plot(v,ce1)
%ce2 = mathieu_ce(2,1,pi/2-v);  % Correct -> is negative of ce1
%hold on; plot(v,ce2) 
%legend('ce1','ce2')



end
