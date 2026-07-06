function aLL=Adversarial_LL(dt,n1,n2,aLL_flag,type_flag)
  % A function that computes the adversarial log-likelihood (aLL), to be 
  % used as a model-differencing metric.
  
  % Output a vector or just the endpoint aLL, depending on type_flag.
  if(strcmpi(type_flag,'sum'))
      sum_fxn=@(x) sum(x);
  elseif(strcmpi(type_flag,'cumsum'))
      sum_fxn=@(x) cumsum(x);
  end
  
  % Integral term.
  N1=sum_fxn(n1)*dt; N1=N1-min(N1);
  N2=sum_fxn(n2)*dt; N2=N2-min(N2);
  
  % Compute aLL, depending on aLL_flag.
  if(strcmpi(aLL_flag,'aLL1'))
      aLL=abs(sum_fxn(abs(log(n1)-log(n2))+abs(N1-N2)));
  elseif(strcmpi(aLL_flag,'aLL2'))
      LL1=sum_fxn(log(n1)-N1);
      LL2=sum_fxn(log(n2)-N2);
      aLL=abs(LL1-LL2);
  end
end