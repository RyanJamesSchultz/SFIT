function [n,N]=EQ_ETIS_Forward(t,Te,Me,Mc,Ti,vi,Vi,params,type_flag)
  % Function that will compute the expected rate (n) and number (N) of 
  % earthquakes as a function of time (t), given the model parameters 
  % (params), injection history (Ti,vi,Vi), and event history (Te,Me).
  % 
  % References:
  % Avouac, J.P., Vrain, M., Kim, T., Smith, J., Ader, T., Ross, Z., & Saarno, T. (2021). A Convolution Model for Earthquake Forecasting Derived from Seismicity Recorded During the ST1 Geothermal Project on Otaniemi Campus, Finland. In Proceedings World Geothermal Congress.
  % Ogata, Y. (1988). Statistical models for earthquake occurrences and residual analysis for point processes. Journal of the American Statistical association, 83(401), 9-27, doi: 10.1080/01621459.1988.10478560.
  % Schultz, R., Ellsworth, W.L., & Beroza, G.C. (2023). An ensemble approach for characterizing trailing-indcued seismicity, Seismological Research Letters, 94(2A), 699–707, doi: 10.1785/0220220352.
  % S&W26. FIX.
  % Shapiro, S. A., Dinske, C., Langenbruch, C., & Wenzel, F. (2010). Seismogenic index and magnitude probability of earthquakes induced during reservoir fluid stimulations. The Leading Edge, 29(3), 304-309, doi: 10.1190/1.3353727.
  % Zhuang, J., Harte, D., Werner, M.J., Hainzl, S., & Zhou, S. (2012). Basic models of seismicity: temporal models, Community Online Resource for Statistical Seismicity Analysis, v1.0, p 42, doi: 10.5078/corssa-79905851.
  % 
  % Written by Ryan Schultz.
  
  % Exponentiate the parameters.
  params([1,3,4,6,8])=10.^params([1,3,4,6,8]);
  
  % Truncate on the magnitude-of-completeness (Mc).
  Te=Te(Me>=Mc);
  Me=Me(Me>=Mc)-Mc;
  
  % Sort the event order based on magnitude, with the smallest first.
  [Me,I]=sort(Me,'ascend');
  Te=Te(I);
  
  % Preallocate space for the output rate & count vectors.
  n=zeros(size(t));
  N=n;
  
  % Compute the flagged seismicity rate/count vectors.
  if(any(strfind(type_flag,'BG')))
      [nb,Nb]=lambda_background(t,Te,Me,Ti,vi,Vi,params); % [Ogata, 1988; Zhuang et al., 2012].
      n=n+nb; N=N+Nb;
  end
  if(any(strfind(type_flag,'ISt')))
      [ni,Ni]=lambda_injections_trail(t,Te,Me,Ti,vi,Vi,params); % [Shapiro et al., 2010; Schultz et al., 2023].
      n=n+ni; N=N+Ni;
  end
  if(any(strfind(type_flag,'ISc')))
      [ni,Ni]=lambda_injections_conv(t,Te,Me,Ti,vi,Vi,params); % [Avouac et al., 2021].
      n=n+ni; N=N+Ni;
  end
  if(any(strfind(type_flag,'ISo')))
      [ni,Ni]=lambda_injections_oracle(t,Te,Me,Ti,vi,Vi,params); % [Schultz & Wiemer, 2026].
      n=n+ni; N=N+Ni;
  end
  if(any(strfind(type_flag,'AS')))
      [na,Na]=lambda_aftershock(t,Te,Me,Ti,vi,Vi,params); % [Ogata, 1988; Zhuang et al., 2012].
      n=n+na; N=N+Na;
  end
  
return
%%% End of main function.




%%% SUBROUTINES.

% Subroutine to compute the background seismicity rate and cumulative counts.
function [n,N]=lambda_background(t,Te,Me,Ti,vi,Vi,params)
  
  % Check if this is required.
  if(any(isnan(params(1))))
      n=0; N=0;
      return;
  end
  
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
  
  % Compute the background seismicity rate and cumulative count.
  n=mu*ones(size(t));
  N=mu*t;
  
return


% Subroutine to compute the injection-induced seismicity rate and cumulative counts (injection & trailing).
function [n,N]=lambda_injections_trail(t,Te,Me,Ti,vi,Vi,params)
  
  % Check if this is required.
  if(any(isnan(params(2:5))))
      n=0; N=0;
      return;
  end
  
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
  %dts=Ti(2)-Ti(1);
  
  % Compute the induced seismicity rate and cumulative event counts (injection vectors).
  ni=vi*10^SI;
  Ni=Vi*10^SI;
  
  % Compute the injection rate and cumulative volume (sampling vectors).
  v=interp1(Ti,vi,t,'nearest',0);
  V=interp1(Ti,Vi,t,'linear','extrap');
  
  % Compute the induced seismicity rate and cumulative event counts (sampling vectors).
  n=v*10^SI;
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
          Te=max([t(end) Ti(end)]);
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
      J2=find((t>=Ts)&(t<=Te));
      if(isempty(J2))
          continue;
      end
      ntj=Kt./((t(J2)+ct-Ts).^pt);
      if(pt==1)
          Ntj=Kt*(log(t(J2)+ct-Ts)-log(ct))+Ns;
      else
          Ntj=(Kt/(pt-1))*(ct^(1-pt)-(t(J2)+ct-Ts).^(1-pt))+Ns;
      end
      n(J2)=ntj;
      N(J2)=Ntj;
      ij=J2(end)+1;
      if(i~=length(I))
          N(ij:end)=N(ij:end)+Nte-Nie;
      end
      
  end
  
return


% Subroutine to compute the injection-induced seismicity rate and cumulative counts (convolutional).
function [n,N]=lambda_injections_conv(t,Te,Me,Ti,vi,Vi,params)
  
  % Check if this is required.
  if(any(isnan(params(2:5))))
      n=0; N=0;
      return;
  end

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
  % Omori-like kernel
  Nk=(ct^(-pt)); % events/day.
  fk=@(tc) ((+1*pt).*((tc+ct).^(-pt-1)))/Nk; % events/day².
  %Fk=@(tc) ((-1   ).*((tc+ct).^(-pt  )))/Nk; % events/day.
  
  % Do the convolution integral numerically.
  nc=(10^SI)*conv(vi,fk(Tk),'full')*dts;
  Tc=linspace(Ti(1),Ti(end)+Tk(end),length(nc));
  Nc=cumsum(nc(1:end-1).*diff(Tc));
  %Nc=conv(v,Nc,'full');

  % Compute the injection rate and cumulative volume (sampling vectors).
  n=interp1(Tc,nc,t,'nearest',0);
  N=interp1(Tc(1:end-1),Nc,t,'linear','extrap');

return


% Subroutine to compute the injection-induced seismicity rate and cumulative counts (ORACLE-like).
function [n,N]=lambda_injections_oracle(t,Te,Me,Ti,vi,Vi,params)
  
  % Check if this is required.
  if(any(isnan(params(2:5))))
      n=0; N=0;
      return;
  end
  
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
  v=interp1(Ti,vi,t,'nearest',0);
  V=interp1(Ti,Vi,t,'linear','extrap');
  
  % Compute the induced seismicity rate and cumulative event counts (sampling vectors).
  n=v*10^SI;
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
          Te=max([t(end) Ti(end)]);
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
      J2=find((t>=Ts)&(t<=Te));
      if(isempty(J2))
          continue;
      end
      ntj=Kt./((t(J2)+ct-Ts).^pt);
      if(pt==1)
          Ntj=Kt*(log(t(J2)+ct-Ts)-log(ct))+Ns;
      else
          Ntj=(Kt/(pt-1))*(ct^(1-pt)-(t(J2)+ct-Ts).^(1-pt))+Ns;
      end
      n(J2)=ntj;
      N(J2)=Ntj;
      ij=J2(end)+1;
      if(i~=length(I))
          N(ij:end)=N(ij:end)+Nte-Nie;
      end
      
  end
  
return


% Subroutine to compute the ETAS seismicity aftershock rate and cumulative counts.
function [n,N]=lambda_aftershock(t,Te,Me,Ti,vi,Vi,params)
  
  % Check if this is required.
  if(any(isnan(params(6:9))))
      n=0; N=0;
      return;
  end
  
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
  
  % Initialize the aftershock rates/counts vectors.
  n=zeros(size(t));
  N=n;
  
  % Loop over the event history.
  for i=1:length(Te)
      
      % (Re)-Initialize the i-th aftershock rates/counts vectors.
      na_i=zeros(size(t));
      Na_i=na_i;
      
      % Get the indicies of relevance, for this event.
      I=(t>Te(i));
      
      % Compute the i-th aftershock rate.
      na_i(I)=exp(alpha*Me(i))*Ka./((t(I)+ca-Te(i)).^pa);
      
      % Compute the i-th aftershock cumulative count.
      if(pa==1)
          Na_i(I)=exp(alpha*Me(i))*Ka*(log(t(I)+ca-Te(i))-log(ca));
      else
          Na_i(I)=exp(alpha*Me(i))*(Ka/(pa-1))*(ca^(1-pa)-(t(I)+ca-Te(i)).^(1-pa));
      end
      
      % Add this i-th aftershock rate/count to the total time series.
      n=n+na_i;
      N=N+Na_i;
  end
  
return


