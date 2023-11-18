function [acov,acor,lag,dof] = autocov(t,y,max_lag)

%   Computes Autocovariance of a time-series
%   after removing linear trend and accounting 
%   for missing data values
%
%   t is time, and assumed to be full. max_lag is the 
%   the maximum lag to compute, in same units as t (i.e., years, months, days, etc)
%   Missing data values must be defined by %nan.
%
%   This version does not interpolate to get full matrices, but
%   uses loops to edit out missing data
%
%   Bin size is determined by the time difference between first 2 points.

  nmax = size(t,1);
  del = t(2) - t(1);
  
  nlags = round(max_lag/del);
  lag = 0:nlags-1;
  lag = lag*del;
  lag = lag';
  epoch = t(1);
  dt = t - epoch;
  
  % check for %nan values 
  
  a = isnan(y);
  nval = 0;
  j = 1;
  for i = 1:nmax
    if (~a(i))
      nval = nval + 1;
      H(j,1) = 1.0;
      H(j,2) = dt(i);
      yt(j,1) = y(i);
      tt(j,1) = dt(i);
      j = j+1;
    end
  end   
      
  nmax2 = size(yt,1);
  
  % Next, estimate and remove linear trend
   
   x = polyfit(tt,yt,1);     
   yfilt = x(2) + x(1)*tt;
   yr = yt - yfilt;
 
    
  acov = zeros(nlags,1);
  ncov = zeros(nlags,1);
  nl = nlags-1;
  
  for i = 1:nmax2
   for j = i-nl:i+nl
     if (j <= nmax2 && j >= 1)
       k = round(abs((tt(j)-tt(i))/del))+1;
       if (k >=1 && k <= nlags)
       acov(k,1) = acov(k,1) + yr(i)*yr(j);
       ncov(k,1) = ncov(k,1) + 1.0;
       end
     end
   end
 end
 
 acov(1) = var(yr);
 for i = 2:nlags
   acov(i) = acov(i)/(ncov(i)-1);
 end

 acor = acov/max(acov);
   
 % compute effective degrees of freedom based on 2x time for acor to first drop 
 % below 1/e
 
 test = 0.0;
 efold = 1/exp(1);
 decor = lag(5);
 for i=2:nlags
     if acor(i) <= efold && test ~= 123.45
         decor = lag(i)*2.0;
         test = 123.45;
     end
 end
 
 ncor = round(decor/del);
 dof = round(nmax2/ncor);
 