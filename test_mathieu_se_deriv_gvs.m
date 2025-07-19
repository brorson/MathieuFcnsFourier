function test_mathieu_se_deriv_gvs(varargin)
  % This reads a file of Mathieu se derivative golden values and uses them
  % to test the output of my se deriv impl.
    
  tol = 2e-1;  % Value is high to accomodate n=34.
               % Also use high value since deriv is computed
	       % using finite differences.
  
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
  
  % The remaining cols hold Mathiue se deriv values for
  % m = 1, 2, ...
  leg = {};
  for i=2:4 % size(M,2)
    sed_gold = M(:,i);
    m = i-1;  % 1, 2, 3, ...
    sed_mine = mathieu_se_deriv(m, q, v)';
    plot(v,sed_mine)
    leg = [leg,['my ',num2str(m)]];    
    hold on
    plot(v,sed_gold)
    leg = [leg,['gv ',num2str(m)]];

    diff(:,i-1) = sed_gold - sed_mine;

  end
  legend(leg)
  title('my sed')
  
  for i=1:size(diff,2)
    ndiff = norm(diff(:,i))/size(diff,2);
    % fprintf('Order = %d, relnormdiff = %e\n', i-1, ndiff)
    if (ndiff > tol)
      fprintf('Failure for order = %d, tol = %e, relnormdiff = %e\n', i-1, tol, ndiff)
    end
  end
  
end
