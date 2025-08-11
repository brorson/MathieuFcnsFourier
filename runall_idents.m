function runall_idents()
  % Runs all identity tests
    
  fcns = {'ce','se','mc1','mc2','ms1','ms2'};
  tests = {'test_mathieu_ce_idents',
	   'test_mathieu_se_idents',
	   'test_mathieu_modmc1_idents',
	   'test_mathieu_modmc2_idents',
	   'test_mathieu_modms1_idents',
	   'test_mathieu_modms2_idents'
	   };
  N = length(fcns);
  pass = zeros(N,1);
  fail = zeros(N,1);
  
  for i=1:N
    fprintf('\n\n')
    fprintf('=====================================\n')
    fprintf('=========   Testing fcn %s  ========\n', fcns{i})
    fh = str2func(tests{i});
    [p, f] = fh();
    pass(i) = p;
    fail(i) = f;
  end
  
  fprintf('=======================================\n')
  fprintf('End of run results\n')
  fprintf('%10s,  %7s,  %7s\n','Fcn','pass','fail')
  for i=1:N
    fprintf('%10s,  %7d,  %7d\n',fcns{i},pass(i),fail(i))
  end
  
end
