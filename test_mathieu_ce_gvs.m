function test_mathieu_ce_gvs()
  % This reads a file of Mathieu ce golden values and uses them
  % to test the output of my ce impl.
    
  
  % Read file holding GVs
  M = csvread('mathieu_ce_gvs.csv');

  % Later I will put the q value into the header of the GV file.
  % Right now I use q = 1;
  q = 1;
  
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
    ndiff = norm(diff(:,i));
    fprintf('Order = %d, normdiff = %e\n', i-1, ndiff)
  end
  
  
end
