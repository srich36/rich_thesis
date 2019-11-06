function boundValue = bound(x,lowerLimit, upperLimit)
  %Return x value between lower and upper limits 
  boundValue=min(max(x,lowerLimit),upperLimit);
end