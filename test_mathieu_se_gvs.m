function test_mathieu_se_gvs()
  % This reads a file of Mathieu se golden values and uses them
  % to test the output of my se impl.
    
  % Later I will put the q value into the header of the GV file.
  % Right now I use q = 1;
  q = 10;
  
  % Read file holding GVs
  %M = csvread('mathieu_se_gvs.csv');
  M = csvread('mathieu_se_gvs_q10.csv');  
  %M = csvread('mathieu_se_gvs_q0.1.csv');    

  % The first col holds the v values.
  v = M(:,1)';
  
  % The remaining cols hold Mathiue se values for
  % m = 1, 2, ...
  leg = {};
  for i=2:size(M,2)
    se_gold = M(:,i);
    m = i-1;
    se_mine = mathieu_se(m, q, v)';
    plot(v,se_mine,'-')
    hold on
    %plot(v,se_gold,'.')
    leg = [leg,num2str(m)];

    diff(:,i-1) = se_gold-se_mine;
  end
  legend(leg)
  title('my se')
  
  for i=1:size(diff,2)
    ndiff = norm(diff(:,i))/size(diff,1);
    fprintf('Order = %d, relnormdiff = %e\n', i, ndiff)
  end
  
  
end
