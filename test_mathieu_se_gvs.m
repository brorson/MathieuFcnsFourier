function test_mathieu_se_gvs(varargin)
  % This reads a file of Mathieu se golden values and uses them
  % to test the output of my se impl.
    
  tol = 2e-6;  % Value is high to accomodate n=34

  % Read filename from the command line if it isn't in the calling args.
  if (length(varargin) == 0)
    fid = stdin();
    filename = fscanf(fid,'%s');
  else
    filename = varargin{1};
  end
  
   fprintf('filename = %s\n', filename)

  % Read data out of the file in the most hacky way possible.
  % The first row holds the q value
  fid = fopen(filename,'r');
  q = fscanf(fid,'%f',1);
  fclose(fid);

  fprintf('Golden value test, filename = %s, q = %f\n', filename, q)
  
  % The remaining rows hold the GV data.
  M = csvread(filename,1,0);  
  
  % The first col holds the v values.
  v = M(:,1)';
  
  % The remaining cols hold Mathiue se values for
  % m = 1, 2, ...
  %leg = {};
  for i=2:size(M,2)
    se_gold = M(:,i);
    m = i-1;
    se_mine = mathieu_se(m, q, v);  % mathieu_se returns col vec.
    %plot(v,se_mine,'-')
    %hold on
    %plot(v,se_gold,'.')
    %leg = [leg,num2str(m)];

    diff(:,i-1) = se_gold-se_mine;
  end
  %legend(leg)
  %title('my se')

  % Iterate over cols which correspond to different Mathieu orders.
  for i=1:size(diff,2)
    % Normalize by number of pts in v since it varies when
    % I modify the GVs.
    ndiff = norm(diff(:,i))/length(v);
    %fprintf('Order = %d, relnormdiff = %e\n', i, ndiff)
    if (ndiff > tol)
      fprintf('Failure for order = %d, tol = %e, relnormdiff = %e\n', i-1, tol, ndiff)
      figure(1)
      plot(v,se_mine,'b-')
      hold on
      plot(v,se_gold,'r.')
      legend('mine','gold')
      title('se')
      %figure(2)
      %plot(v,se_mine-se_gold)
      %title('diff')
      %pause()
      close all;
    end

  end
  
  
end
