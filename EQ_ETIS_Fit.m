function [params,NLL]=EQ_ETIS_Fit(Te,Me,Mc,Ti,vi,Vi,bounds,guess,type_flag,parallel_flag)
  % Function that will MLE fit the model parameters (params), given the 
  % event history (Te,Me) and injection history (Ti,vi,Vi) as data.  Also 
  % the user must specify the observational time bounds (bounds) and 
  % optionally a first guess of the parameters (guess).
  % 
  % Note for units: time should be in days and volumes in m³ (i.e., v in m³/day).
  % Note that SI here isn't the seismogenic index, it's the quantity (SIGMA-b*Mc).
  % Note that Ka is specific to the Mc value input here.
  % 
  % References:
  % Avouac, J.P., Vrain, M., Kim, T., Smith, J., Ader, T., Ross, Z., & Saarno, T. (2021). A Convolution Model for Earthquake Forecasting Derived from Seismicity Recorded During the ST1 Geothermal Project on Otaniemi Campus, Finland. In Proceedings World Geothermal Congress.
  % Ogata, Y. (1988). Statistical models for earthquake occurrences and residual analysis for point processes. Journal of the American Statistical association, 83(401), 9-27, doi: 10.1080/01621459.1988.10478560.
  % Schultz, R., Ellsworth, W.L., & Beroza, G.C. (2023). An Ensemble Approach to Characterizing Trailing‐Induced Seismicity, Seismological Research Letters, 94(2A), 699–707, doi: 10.1785/0220220352.
  % S&W26.  FIX.
  % Shapiro, S. A., Dinske, C., Langenbruch, C., & Wenzel, F. (2010). Seismogenic index and magnitude probability of earthquakes induced during reservoir fluid stimulations. The Leading Edge, 29(3), 304-309, doi: 10.1190/1.3353727.
  % Zhuang, J., Harte, D., Werner, M.J., Hainzl, S., & Zhou, S. (2012). Basic models of seismicity: temporal models, Community Online Resource for Statistical Seismicity Analysis, v1.0, p 42, doi: 10.5078/corssa-79905851.
  % 
  % Written by Ryan Schultz.
  
  % Sort the event order based on magnitude, with the smallest first.
  [Me,I]=sort(Me,'ascend');
  Te=Te(I);
  
  % Truncate on the magnitude-of-completeness (Mc).
  Te=Te(Me>=Mc);
  Me=Me(Me>=Mc)-Mc;
  
  % Define parameter guesses, if none were given.
  if(isempty(guess))
      %     [  mu,   SI,  dt,   ct,   pt,   Ka, alpha,   ca,  pa].
      guess=[1e-4, -2.0, 0.0, 1e-2,  1.2, 5e-3,  5e-1, 1e-2, 1.2];
  else
      % Exponentiate the some of the model parameters.
      guess([1,3,4,6,8])=10.^guess([1,3,4,6,8]);
  end
  
  % Define parameter boundaries.
  %    [  mu,  SI,   dt,    ct,   pt,   Ka, alpha,   ca,   pa].
  Pl = [0.00, -10, +0.0,  1e-4, 1e+0, 0.00,  1e-3, 1e-5, 1e+0];
  Ph = [1e+2, +10, +0.01, 5e-1, 5e+0  5e-1,  1e+1, 1e-0, 5e+0];
  
  % Define the negative log-likelihood function as a Matlab function handle.
  nLL = @(params)  -loglikelihood(Te,Me,Ti,vi,Vi,bounds,params,type_flag);
  
  % Search for the optimal parameters using a Maximum Likelihood Estimator (MLE).
  options = optimoptions(@fmincon,'StepTolerance',1e-4,'ConstraintTolerance',1e-5,'UseParallel',parallel_flag);
  params = fmincon(nLL,guess,[],[],[],[],Pl,Ph,[],options); 
  NLL=nLL(params);
  
  % Zero off any parameters that weren't flagged for fitting.
  if(~any(strfind(type_flag,'BG')))
      params(1)=NaN;
  end
  if(~any(strfind(type_flag,'IS')))
      params(2:5)=NaN;
  end
  if(~any(strfind(type_flag,'AS')))
      params(6:9)=NaN;
  end
  
  % Output some of the model parameters in log-space.
  params([1,3,4,6,8])=log10(params([1,3,4,6,8]));
  
return
%%% End of main function.




%%% SUBROUTINES.

% The log-likelihood function.
function LL=loglikelihood(Te,Me,Ti,vi,Vi,bounds,params,type_flag)
  
  % Preallocate space for the output rate & count vectors.
  n=zeros(size(Te));
  N=[0 0];
  
  % Compute the flagged seismicity rate/count vectors.
  if(any(strfind(type_flag,'BG')))
      [nb,Nb]=lambda_background(Te,bounds,Te,Me,Ti,vi,Vi,params); % Background seismicity [Ogata, 1988; Zhuang et al., 2012].
      n=n+nb; N=N+Nb;
  end
  if(any(strfind(type_flag,'ISt')))
      [ni,Ni]=lambda_injections_trail(Te,bounds,Te,Me,Ti,vi,Vi,params); % Injection-induced seismicity [Shapiro et al., 2010; Schultz et al., 2023].
      n=n+ni; N=N+Ni;
  end
  if(any(strfind(type_flag,'ISc')))
      [ni,Ni]=lambda_injections_conv(Te,bounds,Te,Me,Ti,vi,Vi,params); % Injection-induced seismicity [Avouac et al., 2021].
      n=n+ni; N=N+Ni;
  end
  if(any(strfind(type_flag,'ISo')))
      [ni,Ni]=lambda_injections_oracle(Te,bounds,Te,Me,Ti,vi,Vi,params); % Injection-induced seismicity [Schultz & Wiemer, 2026].
      n=n+ni; N=N+Ni;
  end
  if(any(strfind(type_flag,'AS')))
      [na,Na]=lambda_aftershock(Te,bounds,Te,Me,Ti,vi,Vi,params); % Aftershock seismicity [Ogata, 1988; Zhuang et al., 2012].
      n=n+na; N=N+Na;
  end
  
  % Compute the log-likelihood function [Zhuang et al., 2012].
  LL=sum(log(n))-(max(N)-min(N));
  
return


% Subroutine to compute the background seismicity rate and cumulative counts.
function [n,N]=lambda_background(td,tb,Te,Me,Ti,vi,Vi,params)
  
  % Get the model parameters.
  mu=params(1);
  %SI=params(2);
  %dt=params(3);
  %ct=params(4);
  %pt=params(5);
  %Ka=params(6);
  %alpha=params(7);
  %ca=params(8);
  %pa=params(9);
  
  % Compute the background rate and background cumulative count.
  % Note that rates (n) are evaluated on the data times (td) only and
  % cumulative counts (N) are evaluated on the boundaries times (tb) only.
  n=mu*ones(size(td));
  N=mu*tb;
  
return


% Subroutine to compute the injection-induced seismicity rate and cumulative counts (injection & trailing).
function [n,N]=lambda_injections_trail(td,tb,Te,Me,Ti,vi,Vi,params)
  
  % Get the model parameters.
  %mu=params(1);
  SI=params(2);
  dt=params(3);
  ct=params(4);
  pt=params(5);
  %Ka=params(6);
  %alpha=params(7);
  %ca=params(8);
  %pa=params(9);
  f=1.0;
  
  % Shift the injection time-axis by the time delay.
  Ti=Ti+dt;
  
  % Compute the induced seismicity rate and cumulative event counts (injection vectors).
  ni=vi*10^SI;
  Ni=Vi*10^SI;
  
  % Compute the injection rate and cumulative volume (sampling vectors).
  v=interp1(Ti,vi,td,'nearest',0);
  V=interp1(Ti,Vi,tb,'linear','extrap');
  
  % Compute the induced seismicity rate and cumulative event counts (sampling vectors).
  % Note that rates (n) are evaluated on the data times (td) only and
  % cumulative counts (N) are evaluated on the boundaries times (tb) only.
  n=v*10^SI+eps;
  N=V*10^SI;
  
  % Find all of the time indicies where the induced seismicity rate decreases.
  I=find(diff(ni(:))<0);
  
  % Loop over each of those trailing cases (starting at the end).
  for i=length(I):-1:1
      
      % (Re)-Initialize the j-th trailing rates/counts vectors (injection vectors).
      nt_j=zeros(size(Ti));
      Nt_j=nt_j;
      
      % Get the j-th trailing seismicity start time.
      j1=I(i);
      Ts=Ti(j1);
      %Ts=mean([Ti(j1),Ti(j1+1)]);
      
      % Get the j-th trailing seismicity rate (injection vector).
      J=(Ti>=Ts);
      Kt=f*ni(j1)*ct^pt;
      nt_j(J)=Kt./((Ti(J)+ct-Ts).^pt);
      
      % Get the end time for this trailing seismicity sequence (injection vector).
      j2=(find((nt_j<ni)&(Ti>Ts),1,'first'))-1;
      if(isempty(j2))
          Te=max([max(tb) Ti(end)]);
      else
          Te=Ti(j2+1);
          %Te=mean([Ti(j2),Ti(j2+1)]);
      end
      
      % Find the relevant indicies (injection vectors).
      % Ignore this sequence if there is only one point.
      J=find((Ti>=Ts)&(Ti<Te))+1;
      if(length(J)<=1)
          continue;
      end
      J(J>length(Ti))=[];
      
      % Compute the j-th trailing seismicity cumulative count (injection vector).
      Ns=Ni(j1);
      if(pt==1)
          Nt_j(J)=Kt*(log(Ti(J)+ct-Ts)-log(ct))+Ns;
          Nte    =Kt*(log( Te +ct-Ts)-log(ct))+Ns;
      else
          %size(Ti)
          %size(Nt_j)
          %min(j)
          %max(j)
          %J
          Nt_j(J)=(Kt/(pt-1))*(ct^(1-pt)-(Ti(J)+ct-Ts).^(1-pt))+Ns;
          Nte    =(Kt/(pt-1))*(ct^(1-pt)-( Te +ct-Ts).^(1-pt))+Ns;
      end
      J(end)=[];
      
      % Contribute this trailing seismicity rate/count to the aggregate time series (injection vectors).
      ni(J)=nt_j(J);
      Ni(J)=Nt_j(J);
      
      % Fix the integration at the end time (injection vector).
      ij=J(end)+1;
      Nie=Ni(ij);
      if(i~=length(I))
          Ni(ij:end)=Ni(ij:end)+Nte-Nie;
      end
      
      % Now do this all again, but just for the sampling vectors.
      J2d=find((td>=Ts)&(td<Te));
      if(~isempty(J2d))
          ntj=Kt./((td(J2d)+ct-Ts).^pt);
          n(J2d)=ntj;
      end
      
      J2b=find((tb>=Ts)&(tb<=Te));
      if(~isempty(J2b))
          if(pt==1)
              Ntj=Kt*(log(tb(J2b)+ct-Ts)-log(ct))+Ns;
          else
              Ntj=(Kt/(pt-1))*(ct^(1-pt)-(tb(J2b)+ct-Ts).^(1-pt))+Ns;
          end
      
          %J2b
          N(J2b)=Ntj;
          ij=J2b(end)+1;
          if(i~=length(I))
              N(ij:end)=N(ij:end)+Nte-Nie;
          end
      end
      
  end
  
return


% Subroutine to compute the injection-induced seismicity rate and cumulative counts (convolutional).
function [n,N]=lambda_injections_conv(td,tb,Te,Me,Ti,vi,Vi,params)
  
  % Get the model parameters.
  %mu=params(1);
  SI=params(2);
  dt=params(3);
  ct=params(4);
  pt=params(5);
  %Ka=params(6);
  %alpha=params(7);
  %ca=params(8);
  %pa=params(9);
  
  % Shift the injection time-axis by the time delay.
  Ti=Ti+dt;
  
  % Make the convolution kernel.
  dts=Ti(2)-Ti(1);
  Tk=0:dts:200; % FIX.
  % Omori-like kernel.
  Nk=(ct^(-pt)); % events/day.
  fk=@(tc) ((+1*pt).*((tc+ct).^(-pt-1)))/Nk; % events/day².
  %Fk=@(tc) ((-1   ).*((tc+ct).^(-pt  )))/Nk; % events/day.
  
  % Do the convolution integral numerically.
  nc=(10^SI)*conv(vi,fk(Tk),'full')*dts;
  Tc=linspace(Ti(1),Ti(end)+Tk(end),length(nc));
  Nc=cumsum(nc(1:end-1).*diff(Tc));
  %Nc=conv(v,Nc,'full');
  
  % Compute the injection rate and cumulative volume (sampling vectors).
  n=interp1(Tc,nc,td,'nearest',0)+eps;
  N=interp1(Tc(1:end-1),Nc,tb,'linear','extrap');
  
return


% Subroutine to compute the injection-induced seismicity rate and cumulative counts (ORACLE-like).
function [n,N]=lambda_injections_oracle(td,tb,Te,Me,Ti,vi,Vi,params)
  
  % Get the model parameters.
  %mu=params(1);
  SI=params(2);
  dt=params(3);
  ct=params(4);
  pt=params(5);
  %Ka=params(6);
  %alpha=params(7);
  %ca=params(8);
  %pa=params(9);
  f=1.0;
  
  % Shift the injection time-axis by the time delay.
  Ti=Ti+dt;
  
  % Transform the injection rate into an on/off signal.
  dts=Ti(2)-Ti(1);
  vi=sign(vi>0);
  Vi=cumsum(vi)*dts;
  
  % Compute the induced seismicity rate and cumulative event counts (injection vectors).
  ni=vi*10^SI;
  Ni=Vi*10^SI;
  
  % Compute the injection rate and cumulative volume (sampling vectors).
  v=interp1(Ti,vi,td,'nearest',0);
  V=interp1(Ti,Vi,tb,'linear','extrap');
  
  % Compute the induced seismicity rate and cumulative event counts (sampling vectors).
  % Note that rates (n) are evaluated on the data times (td) only and
  % cumulative counts (N) are evaluated on the boundaries times (tb) only.
  n=v*10^SI+eps;
  N=V*10^SI;
  
  % Find all of the time indicies where the induced seismicity rate decreases.
  I=find(diff(ni(:))<0);
  
  % Loop over each of those trailing cases (starting at the end).
  for i=length(I):-1:1
      
      % (Re)-Initialize the j-th trailing rates/counts vectors (injection vectors).
      nt_j=zeros(size(Ti));
      Nt_j=nt_j;
      
      % Get the j-th trailing seismicity start time.
      j1=I(i);
      Ts=Ti(j1);
      %Ts=mean([Ti(j1),Ti(j1+1)]);
      
      % Get the j-th trailing seismicity rate (injection vector).
      J=(Ti>=Ts);
      Kt=f*ni(j1)*ct^pt;
      nt_j(J)=Kt./((Ti(J)+ct-Ts).^pt);
      
      % Get the end time for this trailing seismicity sequence (injection vector).
      j2=(find((nt_j<ni)&(Ti>Ts),1,'first'))-1;
      if(isempty(j2))
          Te=max([max(tb) Ti(end)]);
      else
          Te=Ti(j2+1);
          %Te=mean([Ti(j2),Ti(j2+1)]);
      end
      
      % Find the relevant indicies (injection vectors).
      % Ignore this sequence if there is only one point.
      J=find((Ti>=Ts)&(Ti<Te))+1;
      if(length(J)<=1)
          continue;
      end
      J(J>length(Ti))=[];
      
      % Compute the j-th trailing seismicity cumulative count (injection vector).
      Ns=Ni(j1);
      if(pt==1)
          Nt_j(J)=Kt*(log(Ti(J)+ct-Ts)-log(ct))+Ns;
          Nte    =Kt*(log( Te +ct-Ts)-log(ct))+Ns;
      else
          %size(Ti)
          %size(Nt_j)
          %min(j)
          %max(j)
          %J
          Nt_j(J)=(Kt/(pt-1))*(ct^(1-pt)-(Ti(J)+ct-Ts).^(1-pt))+Ns;
          Nte    =(Kt/(pt-1))*(ct^(1-pt)-( Te +ct-Ts).^(1-pt))+Ns;
      end
      J(end)=[];
      
      % Contribute this trailing seismicity rate/count to the aggregate time series (injection vectors).
      ni(J)=nt_j(J);
      Ni(J)=Nt_j(J);
      
      % Fix the integration at the end time (injection vector).
      ij=J(end)+1;
      Nie=Ni(ij);
      if(i~=length(I))
          Ni(ij:end)=Ni(ij:end)+Nte-Nie;
      end
      
      % Now do this all again, but just for the sampling vectors.
      J2d=find((td>=Ts)&(td<Te));
      if(~isempty(J2d))
          ntj=Kt./((td(J2d)+ct-Ts).^pt);
          n(J2d)=ntj;
      end
      
      J2b=find((tb>=Ts)&(tb<=Te));
      if(~isempty(J2b))
          if(pt==1)
              Ntj=Kt*(log(tb(J2b)+ct-Ts)-log(ct))+Ns;
          else
              Ntj=(Kt/(pt-1))*(ct^(1-pt)-(tb(J2b)+ct-Ts).^(1-pt))+Ns;
          end
      
          %J2b
          N(J2b)=Ntj;
          ij=J2b(end)+1;
          if(i~=length(I))
              N(ij:end)=N(ij:end)+Nte-Nie;
          end
      end
      
  end
  
return


% Subroutine to compute the ETAS seismicity aftershock rate and cumulative counts.
function [n,N]=lambda_aftershock(td,tb,Te,Me,Ti,vi,Vi,params)
  
  % Get the model parameters.
  %mu=params(1);
  %SI=params(2);
  %dt=params(3);
  %ct=params(4);
  %pt=params(5);
  Ka=params(6);
  alpha=params(7);
  ca=params(8);
  pa=params(9);
  
  % Preallocate space for the output rate & count vectors.
  % Note that rates (n) are evaluated on the data times (td) only and
  % cumulative counts (N) are evaluated on the boundaries times (tb) only.
  n=zeros(size(td))+eps;
  N=zeros(size(tb));
  
  % Loop over the event history.
  for i=1:length(Te)
      
      % (Re)-Initialize.
      na=zeros(size(td));
      Na=zeros(size(tb));
      
      % Get the indicies of relevance, for this event.
      Id=(td>Te(i));
      Ib=(tb>Te(i));
      
      % Compute the i-th aftershock rate.
      na(Id)=exp(alpha*Me(i))*Ka./((td(Id)+ca-Te(i)).^pa);
      
      % Compute the i-th aftershock cumulative count.
      if(pa==1)
          Na(Ib)=exp(alpha*Me(i))*Ka*(log(td(Ib)+ca-Te(i))-log(ca));
      else
          Na(Ib)=exp(alpha*Me(i))*(Ka/(pa-1))*(ca^(1-pa)-(tb(Ib)+ca-Te(i)).^(1-pa));
      end
      
      % Add this i-th epidemic rate/count to the total time series.
      n=n+na;
      N=N+Na;
  end
  
return


