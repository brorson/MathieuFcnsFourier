function test_mathieu_ce_gvs()
  % This reads a file of Mathieu ce golden values and uses them
  % to test the output of my ce impl.
    
  % Later I will put the q value into the header of the GV file.
  % Right now I use q = 1;
  q = 10;
  
  % Read file holding GVs
  %M = csvread('mathieu_ce_gvs.csv');
  M = csvread('mathieu_ce_gvs_q10.csv');    
  %M = csvread('mathieu_ce_gvs_q0.1.csv');  
  
  % The first col holds the v values.
  v = M(:,1)';
  
  % The remaining cols hold Mathiue ce values for
  % m = 0, 1, 2, ...
  leg = {};
  for i=2:size(M,2)
    ce_gold = M(:,i);
    m = i-2;
    ce_mine = mathieu_ce(m, q, v)';
    plot(v,ce_mine)
    hold on
    leg = [leg,num2str(m)];

    diff(:,i-1) = ce_gold-ce_mine;
  end
  legend(leg)
  title('my ce')
  
  for i=1:size(diff,2)
    ndiff = norm(diff(:,i))/size(diff,2);
    fprintf('Order = %d, relnormdiff = %e\n', i-1, ndiff)
  end
  
  
end
