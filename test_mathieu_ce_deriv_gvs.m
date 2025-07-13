function test_mathieu_ce_deriv_gvs()
  % This reads a file of Mathieu ce deriv golden values and uses them
  % to test the output of my ce deriv impl.
    
  % Read file holding GVs
  %filename = 'mathieu_ce_gvs_q10.csv';    
  %filename = 'mathieu_ce_gvs_q0.1.csv';  

  filename = 'mathieu_ce_deriv_gvs_q1.csv';
    
  % The first row holds the q value
  fid = fopen(filename,'r');
  q = fscanf(fid,'%f',1);
  fclose(fid);

  fprintf('Golden value test, filename = %s, q = %f\n', filename, q)
  
  % The remaining rows hold the GV data.
  M = csvread(filename,1,0);  
  
  % The first col holds the v values.
  v = M(:,1)';
  
  % The remaining cols hold Mathiue ce values for
  % m = 0, 1, 2, ...
  leg = {};
  for i=2:size(M,2)
    ced_gold = M(:,i);
    m = i-2;
    ced_mine = mathieu_ce_deriv(m, q, v)';
    plot(v,ced_mine)
    leg = [leg,['my ',num2str(m)]];    
    hold on
    plot(v,ced_gold)
    leg = [leg,['gv ',num2str(m)]];

    diff(:,i-1) = ced_gold - ced_mine;
  end
  legend(leg)
  title('my ced')
  
  for i=1:size(diff,2)
    ndiff = norm(diff(:,i))/size(diff,2);
    fprintf('Order = %d, relnormdiff = %e\n', i-1, ndiff)
  end
  
  
end
