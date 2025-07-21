function test_mathieu_ce_deriv_gvs(varargin)
  % This reads a file of Mathieu ce deriv golden values and uses them
  % to test the output of my ce deriv impl.
    
  tol = 2e-1;  % Value is high to accomodate n=34
               % Also, deriv error is higher since it is computed
	       % using finite differences.

  % Read filename from the command line if it isn't in the calling args.
  if (length(varargin) == 0)
    fid = stdin();
    filename = fscanf(fid,'%s');
  else
    filename = varargin{1};
  end
   
  fprintf('filename = %s\n', filename)

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

  % Create table of differences.
  %leg = {};
  for i=2:size(M,2)
    ced_gold = M(:,i);
    m = i-2;
    ced_mine = mathieu_ce_deriv(m, q, v)';
    %plot(v,ced_mine)
    %leg = [leg,['my ',num2str(m)]];    
    %hold on
    %plot(v,ced_gold)
    %leg = [leg,['gv ',num2str(m)]];

    diff(:,i-1) = ced_gold - ced_mine;
  end
  %legend(leg)
  %title('my ced')
  
  % Now look through table for differences which are too large.
  for i=1:size(diff,2)
    ndiff = norm(diff(:,i))/length(v);
    % fprintf('Order = %d, relnormdiff = %e\n', i-1, ndiff)
    if (ndiff > tol)
      fprintf('Failure for order = %d, tol = %e, relnormdiff = %e\n', i-1, tol, ndiff)
    end
  end
  
  
end
