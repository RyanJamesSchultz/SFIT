% Script to show the ambiguity between induced seismicity rate fits for various model types.
% Used to make Figure S2.
clear;

% Predefine some values.
CaseList={'Basel'}; i=[1];
type_flag_s='ISt+AS'; IS_flag_s='ISt';
type_flag_p='ISc+AS'; IS_flag_p='ISc';
type_flag_m='ISo+AS'; IS_flag_m='ISo';
parallel_flag=true;
sample_flag='simple';
b=1.0;
nm=20;

% Get all of the case data.
D=PreProcData(CaseList,i);

% Get requisite catalogue/injection data.
Meq=D.M';
Teq=D.T';
t=D.t';
v=D.v'; % m³/day.
V=D.V'; % m³.

% Apply Mc filter.
Mc=D.Mc;
I=Meq>=Mc;
Meq=Meq(I);
Teq=Teq(I);

% Temporally filter?
te=t(1)+20000000;
I=(Teq>t(1))&(Teq<te);
Meq=Meq(I);
Teq=Teq(I);

% Start time.
Tbounds=[min(t) max(Teq)];
t=days(t-min(Tbounds));
Teq=days(Teq-min(Tbounds));
te=days(te-min(Tbounds));
Tbounds=days(Tbounds-min(Tbounds));
dT=median(diff(t));

% Trim off the useless parts of the injection time-series.
j=min([find(v>0,1,'last')+3,length(t)]);
t=t(1:j);
v=v(1:j);
V=V(1:j);

% Define the initial guess for the model parameters.
%       [log-mu        SI    log-dt    log-ct        pt    log-Ka     alpha    log-ca        pa]
guess_s=[  -Inf   -1.3778   -2.0000   -0.3010    1.4222   -0.3010    0.3823   -0.1343    3.9389]; % Basel.
guess_p=[  -Inf   -1.3215   -2.2968   -0.6469    1.0152   -0.3527    0.0393   -0.1768    3.9123]; % Basel.
guess_m=[  -Inf    1.2177   -2.2304   -0.3032    1.0005   -0.3337    0.0011   -0.2215    4.9996]; % Basel.
%guess_s=[  -Inf    0.3961   -2.0000   -0.3720    4.9999   -5.3209    0.1382   -5.0000    1.8999]; % FORGE-s3
%guess_p=[  -Inf    0.7135   -2.5069   -0.3806    4.9637   -2.6994    0.2493   -5.0000    1.3402]; % FORGE-s3
%guess_m=[  -Inf    3.6910   -2.0000   -1.7372    1.6165   -3.2051    0.4596   -2.0017    2.6351]; % FORGE-s3

% Fit the model parameters via Maximum Likelihood Estimation.
%tic; [params_s,~]=EQ_ETIS_Fit(Teq,Meq,Mc,t,v,V,bounds,guess_s,type_flag_s,parallel_flag); toc;
%tic; [params_p,~]=EQ_ETIS_Fit(Teq,Meq,Mc,t,v,V,bounds,guess_p,type_flag_p,parallel_flag); toc;
%tic; [params_m,~]=EQ_ETIS_Fit(Teq,Meq,Mc,t,v,V,bounds,guess_m,type_flag_m,parallel_flag); toc;

% Basel fit params.
%        [log-mu        SI    log-dt    log-ct        pt    log-Ka     alpha    log-ca        pa]
params_s=[  -Inf   -1.3778   -2.0000   -0.3010    1.4222   -0.3010    0.3826   -0.1343    3.9389];
params_p=[  -Inf   -1.3215   -2.2968   -0.6469    1.0152   -0.3527    0.0393   -0.1768    3.9123];
params_m=[  -Inf    1.2177   -2.2304   -0.3032    1.0005   -0.3337    0.0011   -0.2216    4.9996];
% FORGE-s3 fit params.
%        [log-mu        SI    log-dt    log-ct        pt    log-Ka     alpha    log-ca        pa]
%params_s=[  -Inf    0.3961   -2.0000   -0.3720    4.9999   -5.3209    0.1382   -5.0000    1.8999];
%params_p=[  -Inf    0.7135   -2.5069   -0.3806    4.9637   -2.6994    0.2493   -5.0000    1.3402];
%params_m=[  -Inf    3.6910   -2.0000   -1.7372    1.6165   -3.2051    0.4596   -2.0017    2.6351];

% Compute the estimated rate/counts.
[n_s,N_s]=EQ_ETIS_Forward(Teq,Teq,Meq,Mc,t,v,V,params_s,type_flag_s); nLL_s=sum(log(n_s))+(max(N_s)-min(N_s));
[n_p,N_p]=EQ_ETIS_Forward(Teq,Teq,Meq,Mc,t,v,V,params_p,type_flag_p); nLL_p=sum(log(n_p))+(max(N_p)-min(N_p));
[n_m,N_m]=EQ_ETIS_Forward(Teq,Teq,Meq,Mc,t,v,V,params_m,type_flag_m); nLL_m=sum(log(n_m))+(max(N_m)-min(N_m));

% Compute the estimated rate/counts from just the injection forcing function.
[n_si,N_si]=EQ_ETIS_Forward(Teq,Teq,Meq,Mc,t,v,V,params_s,IS_flag_s);
[n_pi,N_pi]=EQ_ETIS_Forward(Teq,Teq,Meq,Mc,t,v,V,params_p,IS_flag_p);
[n_mi,N_mi]=EQ_ETIS_Forward(Teq,Teq,Meq,Mc,t,v,V,params_m,IS_flag_m);

% Make the Basel-2 injection scenario.
v2=3*v.*linspace(1,0.5,length(v));
v2(v2>0)=fliplr(v2(v2>0));
V2=cumsum(v2)*dT;

% Make the synthetic catalogues.
[Ts,Ms,Is]=EQ_ETIS_Sample(t,v2,V2,params_s,type_flag_s,b,Mc,Tbounds,dT,sample_flag);
[Tp,Mp,Ip]=EQ_ETIS_Sample(t,v2,V2,params_p,type_flag_p,b,Mc,Tbounds,dT,sample_flag);
[Tm,Mm,Im]=EQ_ETIS_Sample(t,v2,V2,params_m,type_flag_m,b,Mc,Tbounds,dT,sample_flag);

% Estimate the effective seismicity rate.
Rs=movmean(1./diff([0,Ts]),nm);
Rp=movmean(1./diff([0,Tp]),nm);
Rm=movmean(1./diff([0,Tm]),nm);

%%% Plot.
GREY=[0.85,0.85,0.85];
ORANGE='#E6B219';
GREEN='#32CD5A';
TEAL='#27BBD8';

% Time series data.
figure(1); clf;
% Seismicity rates.
ax1a=subplot(211);
semilogy(Teq(1:end-1),1./diff(Teq),'or',MarkerFaceColor='r',HandleVisibility='off');  hold on;
semilogy(t,v,'-b',HandleVisibility='off');
semilogy(Teq,n_s,'-',Color=GREEN,DisplayName='Statistics-based');
semilogy(Teq,n_p,'-',Color=TEAL,DisplayName='Physics-based');
semilogy(Teq,n_m,'-',Color=ORANGE,DisplayName='Machine-learning-based');
semilogy(Teq,n_si,':',Color=GREEN,HandleVisibility='off');
semilogy(Teq,n_pi,':',Color=TEAL,HandleVisibility='off');
semilogy(Teq,n_mi,':',Color=ORANGE,HandleVisibility='off');
xlabel('Time (days)'); ylabel('Injection/earthquake rate (events/day)');
legend(Location='southwest');
% Cumulative event counts.
ax1b=subplot(212);
plot(Teq,(1:length(Teq)),'-r',HandleVisibility='off'); hold on;
plot(Teq,N_s,'-',Color=GREEN,DisplayName='Statistics-based');
plot(Teq,N_p,'-',Color=TEAL,DisplayName='Physics-based');
plot(Teq,N_m,'-',Color=ORANGE,DisplayName='Machine-learning-based');
xlabel('Time (days)'); ylabel('Cumulative earthquakes (counts)');
legend(Location='southeast');
% Link up all of the x-axes.
linkaxes([ax1a,ax1b],'x');
xlim([min(Teq) max(Teq)]);

% Synthetic time series data.
figure(52); clf;
% Seismicity rates.
axS1a=subplot(211);
semilogy(t,v2,'-b',HandleVisibility='off'); hold on;
semilogy(Ts,Rs,'-',Color=GREEN,DisplayName='Statistics-based');
semilogy(Tp,Rp,'-',Color=TEAL,DisplayName='Physics-based');
semilogy(Tm,Rm,'-',Color=ORANGE,DisplayName='Machine-learning-based');
xlabel('Time (days)'); ylabel('Injection/earthquake rate (events/day)');
legend(Location='southwest');
% Cumulative event counts.
axS1b=subplot(212);
plot(Ts,1:length(Ms),'-',Color=GREEN,DisplayName='Statistics-based'); hold on;
plot(Tp,1:length(Mp),'-',Color=TEAL,DisplayName='Physics-based');
plot(Tm,1:length(Mm),'-',Color=ORANGE,DisplayName='Machine-learning-based');
xlabel('Time (days)'); ylabel('Cumulative earthquakes (counts)');
legend(Location='southeast');
% Link up all of the x-axes.
linkaxes([axS1a,axS1b],'x');
xlim([min(Teq) max(Teq)]);

% Report values.
[N_s(end) N_p(end) N_m(end)]-length(Teq)
[nLL_s nLL_p nLL_m]