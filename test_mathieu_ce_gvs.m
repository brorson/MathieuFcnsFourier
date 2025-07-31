function test_mathieu_ce_gvs(varargin)
  % This reads a file of Mathieu ce golden values and uses them
  % to test the output of my ce impl.

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
  
  % The remaining cols hold Mathiue ce values for
  % m = 0, 1, 2, ...
  %leg = {};
  for i=2:size(M,2)
    ce_gold = M(:,i);
    m = i-2;
    ce_mine = mathieu_ce(m, q, v);  % mathieu_ce returns col vector
    %plot(v,ce_mine)
    %hold on
    %leg = [leg,num2str(m)];

    diff(:,i-1) = ce_gold-ce_mine;
  end
  %legend(leg)
  %title('my ce')
  
  % Iterate over cols which correspond to different Mathieu orders.
  for i=1:size(diff,2)
    % Normalize by number of pts in v since it varies when
    % I modify the GVs.
    ndiff = norm(diff(:,i))/length(v);
    %fprintf('Order = %d, relnormdiff = %e\n', i-1, ndiff)
    if (ndiff > tol)
      fprintf('Failure for order = %d, tol = %e, relnormdiff = %e\n', i-1, tol, ndiff)
    end
  end
  
  
end
