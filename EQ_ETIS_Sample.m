function [Ts,Ms,Is]=EQ_ETIS_Sample(t,v,V,params,type_flag,b,Mc,Tbounds,dT, sample_flag)
  % Function that will take fitted ETIS parameters and then sample a 
  % synthetic catalogue.
  %
  % References:
  % Zhuang, J., & Touati, S. (2015). Stochastic simulation of earthquake catalogs. Community Online Resource for Statistical Seismicity Analysis, v1.0, p. 29, doi: 10.5078/corssa-43806322.
  %  
  % Written by Ryan Schultz.

  % Modifiy the seismogenic index for the new Mc.
  %sigma=params(2)+b*Mc1;
  %params(2)=sigma-b*Mc2;
  
  % Modify the aftershock productivity for the new Mc.
  %params(6)=params(6)*exp(params(7)*(Mc2-Mc1));
  
  % Sample the synthetic catalogue, based on the user-defined method.
  if(strcmpi(sample_flag,'simple'))
      [Ts,Ms,Is]=Sample_simple(t,v,V,params,type_flag,b,Mc,Tbounds,dT);
  elseif(strcmpi(sample_flag,'Poisson'))
      [Ts,Ms]=Sample_poisson(t,v,V,params,type_flag,b,Mc,Tbounds,dT);
      Is=[];
  end
  
end




%%% SUBROUTINES.

function [Ts,Ms,Is]=Sample_simple(t,v,V,params,type_flag,b,Mc,Tbounds,dT)
  % This approach uses the inverse transform sampling method.  Its 
  % technically incorrect, as it doesn't give a Poisson process for 
  % inter-event times.
  % That said, its simple and fast.

  % Predefine inputs for GR-MFD.
  params_gr=[0,b];
  Mb=[Mc Inf];
  
  % Make a time vector from start to end.
  Tstart=min(Tbounds);
  Tend=max(Tbounds);
  Tt=Tstart:dT:Tend;

  % Predefine the outputs
  Ts=[];
  Ms=[];
  Is=[];

  % Get the rates/counts along the sampling time vector.
  type_flag2=erase(type_flag,'AS');
  [nt,Nt]=EQ_ETIS_Forward(Tt,Ts,Ms,Mc,t,v,V,params,type_flag2);
  
  % Get the CDF.
  CDF=Nt-Nt(1);
  N=poissrnd(floor(CDF(end)));
  CDF=CDF/CDF(end);
  
  % Precondition to ensure interp1 will work here.
  [CDF,I]=unique(CDF);
  Tt2=Tt(I);
  
  % Sample times from the CDF, via inversion.
  rs=sort(rand([1 N]),'ascend');
  Ts=interp1(CDF,Tt2,rs,'linear','extrap');
  
  % Sample magnitudes.
  [Ms,~]=GR_MFD_Rand(Mb,params_gr,[1 N]);
  Is=zeros(size(Ms));
  
  % Check if we need to account for aftershocks too.
  if(~any(strfind(type_flag,'AS')))
      return
  end
  
  % Prep for loop.
  Ts_i=Ts;
  Ms_i=Ms;
  i=1;
  
  % Loop until there are no more aftershocks needed to account for.
  while(true)
      
      % Get the aftershock rates/counts along the sampling time vector.
      [nt,Nt]=EQ_ETIS_Forward(Tt,Ts_i,Ms_i,Mc,t,v,V,params,'AS');
      
      % Get the CDF.
      CDF=Nt-Nt(1);
      N=poissrnd(floor(CDF(end)));
      CDF=CDF/CDF(end);
      %N
      
      % Check if we can break from the loop.
      if(N==0)
          break
      end
      
      % Precondition to ensure interp1 will work here?
      [CDF,I]=unique(CDF);
      Tt2=Tt(I);
      
      % Sample times from the CDF, via inversion.
      rs=sort(rand([1 N]),'ascend');
      Ts_i=interp1(CDF,Tt2,rs,'linear','extrap');
      
      % Sample magnitudes.
      [Ms_i,~]=GR_MFD_Rand(Mb,params_gr,[1 N]);
      Is_i=i*ones(size(Ms));
      
      % Append new events.
      Ms=[Ms, Ms_i];
      Ts=[Ts, Ts_i];
      Is=[Is, Is_i];
      i=i+1;
      
      % Sort chronologically.
      [Ts,I]=sort(Ts);
      Ms=Ms(I);
      Is=Is(I);
      
  end
end


function [Ts,Ms]=Sample_poisson(t,v,V,params,type_flag,b,Mc,Tbounds,dT)
  % Function that will take MLE fitted parameters and generate a synthetic 
  % catalogue of earthquakes consistent with these parameters.  This code
  % uses the inversion method for a Poisson process [Zhuang et al., 2015].
  
  % Predefine inputs for GR-MFD.
  params_gr=[0,b];
  Mb=[Mc Inf];
  
  % Make a time vector from start to end.
  Tstart=min(Tbounds);
  Tend=max(Tbounds);
  Tt=Tstart:dT:Tend;
  i=1;
  
  % Predefine the outputs
  Ts=[];
  Ms=[];
  
  % Get the rates/counts along the sampling time vector.
  [nt,Nt]=EQ_ETIS_Forward(Tt,Ts,Ms,Mc,t,v,V,params,type_flag);
  
  % Loop until we've reached the catalogue end-time.
  while(true)
      i
      
      % Get the adaptive sampling time vector.
      Tt=getTimeVec(Tt,Tstart,Tend,nt,Nt);
      
      % Get the new rates/counts along the sampling time vector.
      [nt,Nt]=EQ_ETIS_Forward(Tt,Ts,Ms,Mc,t,v,V,params,type_flag);
      
      % Make the Poisson process CDF.
      Ni=Nt-Nt(1);
      CDF=1-exp(-Ni);
      
      % Do some stuff to ensure interp1 will work here?
      [CDF,I]=unique(CDF);
      Tt=Tt(I);
      
      % Check to see if we're at the end-point.
      if(isscalar(CDF))
          break
      end
      
      % Randomly sample event times from the (inverted) CDF.
      rs=rand(1);
      %CDF
      %Tt
      Ts(i)=interp1(CDF,Tt,rs,'linear','extrap');
      
      % Sample the event magnitudes.
      [Ms(i),~]=GR_MFD_Rand(Mb,params_gr,1);
      
      % Check to see if we're past the end-point.
      if(Ts(i)>Tend)
          Ts(i)=[];
          Ms(i)=[];
          break
      end
      
      % Prep for next loop and iterate.
      Tstart=Ts(i);
      Tt=Tstart:dT:Tend;
      i=i+1;
      
  end
  
end


function To=getTimeVec(Ti,Tstart,Tend,ni,Ni)
  % Subroutine to get an adaptive time vector.

  % Determine the most conservative sampling needed.
  dT=1/max(ni);
  %Cf=1e-4;
  %dT=Cf*exp(Ni)./ni;
  
  % Make the time vector.
  To=Tstart:dT:Tend;
  
end