function [M,Mmean]=GR_MFD_Rand(Mb,params,Nr)
  % Function that randomly draws earthquake magnitudes from the 
  % Gutenberg-Richter magnitude-frequency distribution (GR-MFD).
  % 
  % Written by Ryan Schultz.
  
  % Get the upper/lower magnitude bounds.
  m1=min(Mb);
  m2=max(Mb);
  
  % Get the parameters.
  a=params(1);
  b=params(2);

  % The PDF & CDF for a bounded GR-MFD.
  f=b*log(10);
  %x=10.^(-b*(m-m1));
  n=1-10.^(-b*(m2-m1));
  %PDF=(f*x)/n; PDF(m>m2)=NaN; PDF(m<m1)=NaN;
  %CDF=(1-x)/n; CDF(m>m2)=NaN; CDF(m<m1)=NaN;
  
  % Get the random CDF values.
  r=rand(Nr);

  % Use an inverted GR-MFD CDF to map the random CDF values into magnitudes.
  M=m1-log10(1-r*n)/b;
  
  % Mean magnitude of the GR-MFD.
  if(isinf(m2))
      Mmean=(m1*f+1)/f;
  else
      Mmean=((m1*f+1)-(1-n)*(m2*f+1))/(f*n);
  end
  
end