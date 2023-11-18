function [yf,yr,x,err,corrmat,r2,n2] = leastsq2(t,y,epoch,nper,periods)

 % Performs least squares fit to a time-series
 % and returns fit (yf), estimated parameters (x), 
 % and the residuals (yr).  
 %
 % Note: t and periods must be in the 
 % same units (i.e., years, meters, temperature, etc)
 % and both t and y must be column 
 % vectors.
 
 % Missing data values must be defined by %nan.
  
  a = isnan(y);
  n = size(y);
  nmax = max(n);
  dt = t - epoch;
  
  % can only evaluate partials at times when there
  % is an observation
  % Loop through data and accumulate H matrix where
  % data exist
  
  
  j = 1;
  for i = 1:nmax
    if (~a(i))
      H(j,1) = 1.0;
      H(j,2) = dt(i);
      yt(j,:) = y(i);
      nparam = 2;
      if (nper >= 1)
         for k = 1:nper     
           freqp = 2*pi/periods(k);
           H(j,nparam+1) = cos(freqp*dt(i));
           H(j,nparam+2) = sin(freqp*dt(i));
           nparam = nparam + 2;
         end
      end
      j = j+1;
    end
  end   
  
  nmax2 = size(yt,1);
  
  HT = H';
  HTy = HT*yt;
  HTH = HT*H;
  HTHinv = inv(HTH);
  
  x = HTHinv*HTy;
    
% compute residuals and uncertainty

   if (nper == 0)      
     yf = x(1) + x(2)*dt;
   end
   
   if (nper >= 1)
      nparam = 2;
      yf = x(1) + x(2)*dt;
      for k = 1:nper
        freqp = 2*pi/periods(k);
        yf = yf + x(nparam+1)*cos(freqp*dt) + x(nparam+2)*sin(freqp*dt);
        nparam = nparam + 2;
      end
   end
   
  for i = 1:nmax
    if (~a(i))
      yr(i,:) = y(i) - yf(i);
    else
      yr(i,:) = NaN;
    end
  end   
  
  sigma = sqrt(diag(HTHinv));
  
  % compute uncertainty (standard error)
  
    for i   = 1: 2+2*nper
       err(i,:)= nanstd(yr)*sigma(i);
    end
  
  % compute correlation matrix

  for i   = 1: 2+2*nper
    for j   = 1: 2+2*nper
       corrmat(i,j)= HTHinv(i,j)/sigma(i)*sigma(j);
    end
  end
  
  % compute variance explained
  
  r2 = 100*(1 - nanstd(yr)^2/nanstd(y)^2); % r-squared
  %r2 = 100*((nanstd(y)-nanstd(yr))/nanstd(y)); % variance reduced
  
  % number of observations (excluding nan values)
  
  n2 = nmax2;
  
  