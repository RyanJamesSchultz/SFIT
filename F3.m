% Script that makes the schematic SFIT diagram.
% Used to make Figure 3.
clear;

% Predefine some values.
dt=1e-2;
f1i=0.05;
f1s=0.20;
v1p=[1.00 0.75 0.50 0.25];
tf2i=1.00;
tf2p=0.25;
v2i=0.75;
v2p=1.00;
type_flag_p='ISc';
type_flag_s='ISt';
type_flag_m='ISo';

% Predefine the ETIS parameters.
%      [log-mu    SI  log-dt  log-ct    pt  log-Ka  alpha  log-ca     pa]
params=[  -Inf  -2.0   -2.00   -0.85  1.80    -Inf    NaN    -Inf    NaN];
Mc=0.0;

% Define the time axis.
t=0:dt:10.0;

% Get params
ct=10^params(4);
pt=params(5);

% Define time scales for SFIT Phase 1.
dT1i=ct*(f1i^(-1/pt)-1);
dT1s=ct*(f1s^(-1/pt)-1);

% Initialize and loop for SFIT Phase 1.
v=zeros(size(t)); t1e=0;
for i=1:length(v1p)
    
    % Iterate some of the Phase 1 pulse details.
    t1s=t1e+2*dT1s;
    t1e=t1s+dT1i;
    
    % Define the injection protocol for one Phase 1 pulse.
    v((t>=t1s)&(t<=t1e))=10^v1p(i);
end

% Define time scales for SFIT Phase 2.
dT2i=dT1i*tf2i;
dT2r=1.1; %FIX?
dT2p=dT1i*tf2p;

% Define the Phase 2 ramp & pulse details.
t2s=6;
t2m=t2s+dT2i;
t2r=t2m+dT2r;
t2e=t2r+dT2p;

% Define the injection protocol for the Phase 2 ramp & pulse.
v((t>=t2s)&(t<=t2m))=10^v2i;
v((t>=t2m)&(t<=t2r))=10.^linspace(v2i, -3, length(v((t>=t2m)&(t<=t2r))) );
v((t>=t2r)&(t<=t2e))=10^v2p;

% Compute the cumulative injected volume.
V=cumsum(10.^v)*dt;

% Compute the three forecast paradigm seismicity rates for SFIT.
[np,Np]=EQ_ETIS_Forward(t,[],[],Mc,t,v,V,params,type_flag_p); np=np*10/max(np);
[ns,Ns]=EQ_ETIS_Forward(t,[],[],Mc,t,v,V,params,type_flag_s); ns=ns*10/max(ns);
[nm,Nm]=EQ_ETIS_Forward(t,[],[],Mc,t,v,V,params,type_flag_m); nm=nm*10/max(nm);




%%% Plot.
ORANGE='#E6B219';
GREEN='#32CD5A';
TEAL='#27BBD8';

%Figure 3.
figure(3); clf;
plot(t,log10(v+1e-6),'-b'); hold on;
plot(t,log10(nm+eps),'-',Color=ORANGE);
plot(t,log10(np+eps),'-',Color=TEAL);
plot(t,log10(ns+eps),'-',Color=GREEN);
plot(xlim(), [-1 -1],':k');
xlabel('Time'); ylabel('Rate');
ylim([-1.25 +1.25]);
xlim([min(t) max(t)]);


