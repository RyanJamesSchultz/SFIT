% Plots the proposed injection protocols from the model-differencing process.
% See also script_MakeDiffModels.m.
% Used to visually examine proposed injection protocols.
clear;

% Define some values.
search_path='Vdata/V_SFIT_sp_aLL1_*.mat';
dt=1e-3;
Ns=100;
plot_type='step';
Ntop=20;

% Basel fit params.                [log-mu        SI    log-dt    log-ct        pt    log-Ka     alpha    log-ca        pa]
param_flag_s='BG+ISt+AS'; params_s=[  -Inf   -1.3778   -2.0000   -0.3010    1.4222   -0.3010    0.3826   -0.1343    3.9389];
param_flag_p='BG+ISc+AS'; params_p=[  -Inf   -1.3215   -2.2968   -0.6469    1.0152   -0.3527    0.0393   -0.1768    3.9123];
param_flag_m='BG+ISo+AS'; params_m=[  -Inf    1.2177   -2.2304   -0.3032    1.0005   -0.3337    0.0011   -0.2216    4.9996];
Mc=0.90;
% FORGE-s3 fit params.
%                                  [log-mu        SI    log-dt    log-ct        pt    log-Ka     alpha    log-ca        pa]
%param_flag_s='BG+ISt+AS'; params_s=[  -Inf    0.3961   -2.0000   -0.3720    4.9999   -5.3209    0.1382   -5.0000    1.8999];
%param_flag_p='BG+ISc+AS'; params_p=[  -Inf    0.7135   -2.5069   -0.3806    4.9637   -2.6994    0.2493   -5.0000    1.3402];
%param_flag_m='BG+ISo+AS'; params_m=[  -Inf    3.6910   -2.0000   -1.7372    1.6165   -3.2051    0.4596   -2.0017    2.6351];
%Mc=-1.20;

% Make the time axis.
ti=0:dt:1-dt;
t=[ti+0,ti+1,ti+2];
it=length(ti)+1;

% Load the data structures and sort based on aLL.
Soln=struct('v',[],'aLL',[]);
input_files=dir(search_path);
for i=1:length(input_files)
    load(fullfile(input_files(i).folder,input_files(i).name),'soln');
    Soln(i)=soln(end);
end
[~,I]=sort([Soln.aLL],'descend'); Soln=Soln(I);
Soln=Soln(1:Ntop);

% Regenerate the model forecasts.
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

% Loop over all the SFIT proposals.
figure(3); clf;
figure(5); clf;
for i=1:length(Soln)
    
    % Two types of aLL.
    aLL1=Adversarial_LL(dt,Soln(i).n1+eps,Soln(i).n2+eps,'aLL1','cumsum'); % Piecewise absolute differences.
    aLL2=Adversarial_LL(dt,Soln(i).n1+eps,Soln(i).n2+eps,'aLL2','cumsum'); % End absolute difference.
    aLL3=(aLL1+aLL2)/2; % Average.

    % Find time to peak.
    %t_peak=t(find((diff([0,aLL2])<0)&(t>1.2),1));
    t_peak=0;

    % Injection protocol proposals.
    figure(3);
    semilogy(t-t_peak,Soln(i).v+1/Ns*10,'-b',HandleVisibility='off'); hold on;
    semilogy(t-t_peak,Soln(i).n1,'-r',DisplayName=Soln(i).model_name1);
    semilogy(t-t_peak,Soln(i).n2,'-m',DisplayName=Soln(i).model_name2);
    title([parts{3},'-',parts{4}]);
    xlabel('Relative Time'); ylabel('Relative Rate');
    
    % Plot the two types of aLL.
    figure(5);
    plot(t-t_peak,aLL1,'-b'); hold on;
    plot(t-t_peak,aLL2,'-c');
    plot(t-t_peak,aLL3,':k');
    title([parts{3},'-',parts{4}]);
    xlabel('Relative Time'); ylabel('Relative Rate');
    if(strcmpi(plot_type,'step'))
        i
        pause
        figure(3); clf;
        figure(5); clf;
    end
end





