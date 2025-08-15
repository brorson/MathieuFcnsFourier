function plot_eigs()
  % Make plot of Mathieu eigs vs. q to reproduce plot
  % on https://dlmf.nist.gov/28.2
  
  % Number of sample points
  N = 1000;
  
  % Domain of q values to examine (for plotting)
  qs = linspace(0,10,N)';
  
  % Number of each type of eigenvalue to track
  Ne = 6;  % Ne of a and Ne of b


  %---------------------------------------------
  % First plot ce eigs

  % Preallocate a vector to store values.
  as = zeros(length(qs), Ne);

  fprintf('Calculating a eigenvalues ... \n')

  % Loop over order m
  for m=0:Ne-1
    % Loop over qs.
    for i = 1:length(qs)
      q = qs(i);  % Get this value of q.
      as(i,m+1) = mathieu_a(m, q);
    end
  end
  
  %fprintf('Even eigs close to q=0:\n')
  %disp(as(1,:)')

  % Make plot of eigenvalues vs. q
  figure(1)
  c = {};
  for j=1:Ne
    hold on
    if (j==1)
      h1 = plot(qs,as(:,j),'b-','LineWidth',2);
      c = [c, ['a',num2str(j-1)]];
    else
      c = [c, ''];
      plot(qs,as(:,j),'b-','LineWidth',2);
    end
  end


  %---------------------------------------------
  % Next plot se eigs
  
  % Preallocate b vector to store values.
  bs = zeros(length(qs), Ne);
  
  fprintf('Calculating b eigenvalues ... \n')

  % Loop over order m
  for m=1:Ne
    % Loop over qs.
    for i = 1:length(qs)
      q = qs(i);  % Get this value of q.
      bs(i,m) = mathieu_b(m, q);
    end
  end
  
  %fprintf('Odd eigs close to q=0:\n')
  %disp(bs(1,:)')

  % Make plot of eigenvalues vs. q
  %figure(1)
  for j=1:Ne
    hold on
    if (j==1)
      c = [c, ['b',num2str(j)]];
      h2 = plot(qs,bs(:,j),'r--','LineWidth',2);
    else
      c = [c, ''];
      plot(qs,bs(:,j),'r--','LineWidth',2);
    end
  end
  
  % Turn these on to reproduce the DLMF plot
  xlim([0,10]);
  ylim([-5,20]);
  

  title('First Mathieu eigenvalues vs. q')
  xlabel('q')
  ylabel('eigenvalue')
  %legend([1,Ne+1],'ce eigs','se eigs')
  %ylim([-50,90])
  %ylim([-100,300])  
  legend([h1,h2],'ce eigs','se eigs','Location','NorthWest')
  
end