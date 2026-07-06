% Plots the proposed injection protocols from the model-differencing process.
% See also script_MakeDiffModels.m.
% Used to make Figures 2 & Sx-SX.
clear;

% Predefine some values.
search_path_sp='Vdata/V_SFIT_sp_aLL1_*.mat';
search_path_sm='Vdata/V_SFIT_sm_aLL1_*.mat';
search_path_pm='Vdata/V_SFIT_pm_aLL1_*.mat';
dt=1e-3;
Ns=100;
plot_type='all';
Ntop=15;

% Make the time axis.
ti=0:dt:1-dt;
t=[ti+0,ti+1,ti+2];
it=length(ti)+1;

% Initialize the data structures.
Soln_sp=struct('v',[],'aLL',[]);
Soln_sm=struct('v',[],'aLL',[]);
Soln_pm=struct('v',[],'aLL',[]);

% Load the data structures.
input_files_sp=dir(search_path_sp);
for i=1:length(input_files_sp)
    load(fullfile(input_files_sp(i).folder,input_files_sp(i).name),'soln');
    Soln_sp(i)=soln(end);
end
input_files_sm=dir(search_path_sm);
for i=1:length(input_files_sm)
    load(fullfile(input_files_sm(i).folder,input_files_sm(i).name),'soln');
    Soln_sm(i)=soln(end);
end
input_files_pm=dir(search_path_pm);
for i=1:length(input_files_pm)
    load(fullfile(input_files_pm(i).folder,input_files_pm(i).name),'soln');
    Soln_pm(i)=soln(end);
end

% Sort the data structures based on aLL.
[~,I]=sort([Soln_sp.aLL],'descend'); Soln_sp=Soln_sp(I);
[~,I]=sort([Soln_sm.aLL],'descend'); Soln_sm=Soln_sm(I);
[~,I]=sort([Soln_pm.aLL],'descend'); Soln_pm=Soln_pm(I);
Soln_sp=Soln_sp(1:Ntop);
Soln_sm=Soln_sm(1:Ntop);
Soln_pm=Soln_pm(1:Ntop);

% Regenerate the model forecasts.
Soln_sp=RegenForecasts(Soln_sp,input_files_sp,t,dt);
Soln_sm=RegenForecasts(Soln_sm,input_files_sm,t,dt);
Soln_pm=RegenForecasts(Soln_pm,input_files_pm,t,dt);


%%% Plot.

% Figure 2.
figure(2); clf;

% Panel a).
subplot(311);
for i=1:length(Soln_sp)
    
    % Two types of aLL.
    aLL1=Adversarial_LL(dt,Soln_sp(i).n1+eps,Soln_sp(i).n2+eps,'aLL1','cumsum'); % Piecewise absolute differences.
    aLL2=Adversarial_LL(dt,Soln_sp(i).n1+eps,Soln_sp(i).n2+eps,'aLL2','cumsum'); % End absolute difference.
    aLL3=(aLL1+aLL2)/2; % Average.

    % Find time to peak.
    %t_peak=t(find((diff([0,aLL2])<0)&(t>1.2),1))-1;
    t_peak=1;

    % Injection protocol proposals.
    v=Soln_sp(i).v+1/Ns*10; v(Soln_sp(i).v==0)=1e-2;
    semilogy(t-t_peak,v,'-b',HandleVisibility='off'); hold on;
    semilogy(t-t_peak,Soln_sp(i).n1,'-r',DisplayName=Soln_sp(i).model_name1);
    semilogy(t-t_peak,Soln_sp(i).n2,'-m',DisplayName=Soln_sp(i).model_name2);
end
title('SP');
xlabel('Relative Time'); ylabel('Relative Rate');
xlim([-0.5 2.5]);

% Panel b).
subplot(312);
for i=1:length(Soln_sm)
    
    % Two types of aLL.
    aLL1=Adversarial_LL(dt,Soln_sm(i).n1+eps,Soln_sm(i).n2+eps,'aLL1','cumsum'); % Piecewise absolute differences.
    aLL2=Adversarial_LL(dt,Soln_sm(i).n1+eps,Soln_sm(i).n2+eps,'aLL2','cumsum'); % End absolute difference.
    aLL3=(aLL1+aLL2)/2; % Average.

    % Find time to peak.
    t_peak=t(find((diff([0,aLL2])<0)&(t>1.2),1))-1;
    %t_peak=1;

    % Injection protocol proposals.
    v=Soln_sm(i).v+1/Ns*10; v(Soln_sm(i).v==0)=1e-2;
    semilogy(t-t_peak,v,'-b',HandleVisibility='off'); hold on;
    semilogy(t-t_peak,Soln_sm(i).n1,'-r',DisplayName=Soln_sm(i).model_name1);
    semilogy(t-t_peak,Soln_sm(i).n2,'-m',DisplayName=Soln_sm(i).model_name2);
end
title('SM');
xlabel('Relative Time'); ylabel('Relative Rate');
xlim([-0.5 2.5]);

% Panel c).
subplot(313);
for i=1:length(Soln_pm)
    
    % Two types of aLL.
    aLL1=Adversarial_LL(dt,Soln_pm(i).n1+eps,Soln_pm(i).n2+eps,'aLL1','cumsum'); % Piecewise absolute differences.
    aLL2=Adversarial_LL(dt,Soln_pm(i).n1+eps,Soln_pm(i).n2+eps,'aLL2','cumsum'); % End absolute difference.
    aLL3=(aLL1+aLL2)/2; % Average.

    % Find time to peak.
    t_peak=t(find((diff([0,aLL2])<0)&(t>1.2),1))-0.5;
    %t_peak=1;

    % Injection protocol proposals.
    v=Soln_pm(i).v+1/Ns*10; v(Soln_pm(i).v==0)=1e-2;
    semilogy(t-t_peak,v,'-b',HandleVisibility='off'); hold on;
    semilogy(t-t_peak,Soln_pm(i).n1,'-r',DisplayName=Soln_pm(i).model_name1);
    semilogy(t-t_peak,Soln_pm(i).n2,'-m',DisplayName=Soln_pm(i).model_name2);
end
title('PM');
xlabel('Relative Time'); ylabel('Relative Rate');
xlim([-0.5 2.5]);




%%%% SUBROUNTINES.

% Function to regenerate the forecasts and stuff the results into the data structure.
function [Soln]=RegenForecasts(Soln,input_files,t,dt)
  
  % Basel fit params.                [log-mu        SI    log-dt    log-ct        pt    log-Ka     alpha    log-ca        pa]
  param_flag_s='BG+ISt+AS'; params_s=[  -Inf   -1.3778   -2.0000   -0.3010    1.4222   -0.3010    0.3826   -0.1343    3.9389];
  param_flag_p='BG+ISc+AS'; params_p=[  -Inf   -1.3215   -2.2968   -0.6469    1.0152   -0.3527    0.0393   -0.1768    3.9123];
  param_flag_m='BG+ISo+AS'; params_m=[  -Inf    1.2177   -2.2304   -0.3032    1.0005   -0.3337    0.0011   -0.2216    4.9996];
  Mc=0.90;
  % FORGE-s3 fit params.
  %        [log-mu        SI    log-dt    log-ct        pt    log-Ka     alpha    log-ca        pa]
  %params_s=[  -Inf    0.3961   -2.0000   -0.3720    4.9999   -5.3209    0.1382   -5.0000    1.8999];
  %params_p=[  -Inf    0.7135   -2.5069   -0.3806    4.9637   -2.6994    0.2493   -5.0000    1.3402];
  %params_m=[  -Inf    3.6910   -2.0000   -1.7372    1.6165   -3.2051    0.4596   -2.0017    2.6351];
  %Mc=-1.20;
  
  % Loop over 
  for i=1:length(Soln)
      parts=split(input_files(i).name,'_');
      model_flag1=parts{3}(1);
      model_flag2=parts{3}(2);
      Soln(i).V=cumsum(Soln(i).v)*dt;
      if(strcmpi(model_flag1,'s'))
          [n1,N1]=EQ_ETIS_Forward(t,[],[],Mc,t,Soln(i).v,Soln(i).V,params_s,param_flag_s);
          Soln(i).model_name1='Statistics-based';
      elseif(strcmpi(model_flag1,'p'))
          [n1,N1]=EQ_ETIS_Forward(t,[],[],Mc,t,Soln(i).v,Soln(i).V,params_p,param_flag_p);
          Soln(i).model_name1='Physics-based';
      elseif(strcmpi(model_flag1,'m'))
          [n1,N1]=EQ_ETIS_Forward(t,[],[],Mc,t,Soln(i).v,Soln(i).V,params_m,param_flag_m);
          Soln(i).model_name1='Machine-learning-based';
      end
      if(strcmpi(model_flag2,'s'))
          [n2,N2]=EQ_ETIS_Forward(t,[],[],Mc,t,Soln(i).v,Soln(i).V,params_s,param_flag_s);
          Soln(i).model_name2='Statistics-based';
      elseif(strcmpi(model_flag2,'p'))
          [n2,N2]=EQ_ETIS_Forward(t,[],[],Mc,t,Soln(i).v,Soln(i).V,params_p,param_flag_p);
          Soln(i).model_name2='Physics-based';
      elseif(strcmpi(model_flag2,'m'))
          [n2,N2]=EQ_ETIS_Forward(t,[],[],Mc,t,Soln(i).v,Soln(i).V,params_m,param_flag_m);
          Soln(i).model_name2='Machine-learning-based';
      end
      Soln(i).n1=n1/max(N1); Soln(i).N1=N1/max(N1);
      Soln(i).n2=n2/max(N2); Soln(i).N2=N2/max(N2);
  end
end