clear;

% Predefine some values.
CaseList={'FORGE22'}; i=[3];
BG_flag='';
IS_flag='ISt';
AS_flag='AS';
parallel_flag=true;
dTs=0.001;

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
bounds=[min(t) max(Teq)];
t=days(t-min(bounds));
Teq=days(Teq-min(bounds));
te=days(te-min(bounds));
bounds=days(bounds-min(bounds));

% Trim off the useless parts of the injection time-series.
j=min([find(v>0,1,'last')+3,length(t)]);
t=t(1:j);
v=v(1:j);
V=V(1:j);

% Define the initial guess for the model parameters.
%      log-mu,     SI, log-dt, log-ct,   pt, log-Ka,  alpha, log-ca,   pa].
guess=[-4.000, -2.400, -1.000, -0.500, 1.50, -1.100, 0.3500, -0.500, 3.10];
type_flag=[BG_flag,'+',IS_flag,'+',AS_flag];

%%
% Get the MLE fitted model parameters.
tic;
[params,nLL]=EQ_ETIS_Fit(Teq,Meq,Mc,t,v,V,bounds,guess,type_flag,parallel_flag);
toc;

% Compute the estimated rate/counts.
[n1,N1]=EQ_ETIS_Forward(Teq,Teq,Meq,Mc,t,v,V,params,type_flag);
[nb,Nb]=EQ_ETIS_Forward(Teq,Teq,Meq,Mc,t,v,V,params,BG_flag);
[ni,Ni]=EQ_ETIS_Forward(Teq,Teq,Meq,Mc,t,v,V,params,IS_flag);
[na,Na]=EQ_ETIS_Forward(Teq,Teq,Meq,Mc,t,v,V,params,AS_flag);

% Compute CDF residuals.
r1=residuals1(Teq,Teq,N1);

% Compute the running log-likelihoods.
nll1=(cumsum(log(n1))-(N1-N1(1)));

%%
%%% Plot.
GREY=[0.85,0.85,0.85];

% Time series data.
figure(1); clf;
% Raw rates.
ax1a=subplot(411);
semilogy(Teq(1:end-1),1./diff(Teq),'or');  hold on;
semilogy(t,v,'-b');
semilogy(Teq,n1,'-m');
xlabel('Time (days)'); ylabel('Injection/earthquake rate (events/day)');
% Binned rates.
ax1b=subplot(412);
histogram(Teq,min(bounds):dTs:max(bounds),FaceColor='r'); hold on;
plot(Teq,n1*dTs,'-m');
xlabel('Time (days)'); ylabel('Binned earthquake rate (events/bin)');
set(gca, 'YScale', 'log');
% Cumulative counts.
ax1c=subplot(413);
plot(Teq,N1,'-m'); hold on;
plot(Teq,(1:length(Teq)),'-r');
xlabel('Time (days)'); ylabel('Cumulative earthquakes (counts)');
% MvT.
ax1d=subplot(414);
plot(Teq,Meq,'or');
xlabel('Time (days)'); ylabel('Magnitude');
linkaxes([ax1a,ax1b,ax1c,ax1d],'x');
xlim([min(Teq) max(Teq)]);

% Plot CDF residuals and GoF metrics.
figure(2); clf;
% Running nLL (linear-scale)
ax2a=subplot(211);
plot(Teq,nll1./(1:length(nll1)),'-m'); hold on;
plot(Teq(2:end),diff(nll1),'dm');
plot(te*[1 1], ylim(),'--k');
xlabel('Time (days)'); ylabel('Running (normalized) log-likelihood');
ylim([-5 10]);
% Residual CDFs (linear-scale)
Ylim=0.5*[-1 +1];
ax2b=subplot(212);
plot(Teq,zeros(size(Teq)),'-r'); hold on;
plot(Teq,r1,'-m');
plot(te*[1 1], Ylim,'--k');
xlabel('Time (days)'); ylabel('Cumulative Distribution Residual (-)');
ylim(Ylim);

% Link up all of the x-axes.
linkaxes([ax1a,ax1b,ax1c,ax1d,ax2a,ax2b],'x');
xlim([min(Teq) max(Teq)]);

% Report values.
params
length(Teq)-N1(end)
[Nb(end) Ni(end) Na(end)]/N1(end)
nLL
nLL/length(Teq)




% SUBROUTINES.

% Compute the CDF residual.
function res = residuals1(xd,xf,yf)
  % Calculates the CDF residual (on fitted x-axis).
  yd=interp1(xd+1e-10*(1:length(xd)),1:length(xd),xf,'linear','extrap');
  res=(yf-yd)/length(xd);
end