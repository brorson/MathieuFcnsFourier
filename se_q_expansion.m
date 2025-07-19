function ce = se_q_expansion(m,q,v)
  % This implements the general expansion for small q shown
  % in DLMF 28.6.26
    
  t1 = sin(m*v);
  t2 = q*(sin((m+2)*v)/(m+1) - sin((m-2)*v)/(m-1))/4;
  t3 = (q^2)*(sin((m+4)*v)/((m+1)*(m+2)) + sin((m-4)*v)/((m-1)*(m-2)) - 2*(m*m+1)/((m*m-1)^2)*sin(m*v))/32;
  ce = t1-t2+t3;
end