% A script that will use Markov chain Monte Carlo (McMC) to maximize the
% adversarial log-likelihood (aLL) between to injection forcing functions.
clear;

% Define some values.
dt=1e-3;
Ns=100;
model_flag1='p';
model_flag2='m';
aLL_flag='aLL1';
Niter=1000;
Nr=1;

% Basel fit params.                [log-mu        SI    log-dt    log-ct        pt    log-Ka     alpha    log-ca        pa]
param_flag_s='BG+ISt+AS'; params_s=[  -Inf   -1.3778   -2.0000   -0.3010    1.4222   -0.3010    0.3826   -0.1343    3.9389];
param_flag_p='BG+ISc+AS'; params_p=[  -Inf   -1.3214   -2.2968   -0.6469    1.0152   -0.3527    0.0393   -0.1768    3.9123];
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

% Loop over the McMC realizations.
for l=1:Nr
    
    % Initialize.
    soln=struct('v',[],'aLL',[]);
    j=1; I0=0; L0=0;
    
    % Loop over the McMC iterations.
    for i=1:Niter
        [l, i, j]
        
        % Invent an injection protocol.
        [vi,I1]=RandomInjection(ti,Ns,I0);
        v1=[zeros(size(ti)),vi,zeros(size(ti))];
        V1=cumsum(v1)*dt;
        
        % Forecast the seismicity response for each model (and normalize).
        if(strcmpi(model_flag1,'s'))
            [n1,N1]=EQ_ETIS_Forward(t,[],[],Mc,t,v1,V1,params_s,param_flag_s);
            model_name1='Statistics-based';
        elseif(strcmpi(model_flag1,'p'))
            [n1,N1]=EQ_ETIS_Forward(t,[],[],Mc,t,v1,V1,params_p,param_flag_p);
            model_name1='Physics-based';
        elseif(strcmpi(model_flag1,'m'))
            [n1,N1]=EQ_ETIS_Forward(t,[],[],Mc,t,v1,V1,params_m,param_flag_m);
            model_name1='Machine-learning-based';
        end
        if(strcmpi(model_flag2,'s'))
            [n2,N2]=EQ_ETIS_Forward(t,[],[],Mc,t,v1,V1,params_s,param_flag_s);
            model_name2='Statistics-based';
        elseif(strcmpi(model_flag2,'p'))
            [n2,N2]=EQ_ETIS_Forward(t,[],[],Mc,t,v1,V1,params_p,param_flag_p);
            model_name2='Physics-based';
        elseif(strcmpi(model_flag2,'m'))
            [n2,N2]=EQ_ETIS_Forward(t,[],[],Mc,t,v1,V1,params_m,param_flag_m);
            model_name2='Machine-learning-based';
        end
        n1=n1/max(N1); N1=N1/max(N1);
        n2=n2/max(N2); N2=N2/max(N2);
        
        % Compute the adversarial log-likelihood metric.
        aLL=Adversarial_LL(dt,n1+eps,n2+eps,aLL_flag,'sum');
        
        % Check if the new injection protocol is worse.
        if( aLL<L0 )
            % If it's worse, roll a saving throw.
            if( (aLL-L0)<log(rand()) )
                continue; % Reject the protocol, if the saving throw is failed.
            end
        end
        
        % Save the injection protocol and iterate.
        soln(j).v=v1; soln(j).aLL=aLL;
        j=j+1; v0=v1; L0=aLL; I0=I1;
    end
    
    % Save the output data structure.
    output_filename=['V_SFIT_',model_flag1,model_flag2,'_',aLL_flag,'_',num2str(l,'%02.f'),'.mat'];
    save(output_filename,'soln');
end


%%%
% Plot.
if(Nr==1)
    
    % Proposed injection protocol.
    figure(3); clf;
    semilogy(t,v1+1/Ns*10,'-b',HandleVisibility='off'); hold on;
    semilogy(t,n1,'-r',DisplayName=model_name1);
    semilogy(t,n2,'-m',DisplayName=model_name2);
    xlabel('Time'); ylabel('Rate');
    legend();
    
    % Adversarial log-likelihood versus McMC iteration.
    figure(4); clf;
    plot(1:length(soln),[soln.aLL],'-xb');
    xlabel('Injection Protocol Iteration'); ylabel('Adversarial Log-Likelihood');
    
    % Two types of aLL.
    aLL1=Adversarial_LL(dt,n1+eps,n2+eps,'aLL1','cumsum'); % Piecewise absolute differences.
    aLL2=Adversarial_LL(dt,n1+eps,n2+eps,'aLL2','cumsum'); % End absolute difference.
    aLL3=(aLL1+aLL2)/2; % Average.
    
    % Plot the two types of aLL.
    figure(5); clf;
    plot(t,aLL1,'-b'); hold on;
    plot(t,aLL2,'-c');
    plot(t,aLL3,':k');
    xlabel('Time'); ylabel('Rate');
    
end






%%%% SUBROUNTINES.

% Invent an injection protocol.
function [v,I]=RandomInjection(t,Ns,I)
  
  % Predefine some parameters.
  Nt=length(t);
  
  % Get/revise a vector of random integers.
  if(I==0)
      I=randi(Nt,[1 Nt*Ns]);
  else
      I=revise(I,Nt,Ns);
  end
  
  % Make the injection protocol.
  v=histcounts(I,Nt)/Ns;
  
end


% Revise a vector of random integers.
function [I]=revise(I,Nt,Ns)
  
  % Randomly decide the remove/add parameters.
  c=2;
  Nr=ceil(abs(normrnd(0.10,0.05))*Nt*Ns);
  Di=randi(c);
  Dl=randi(Nt,1);
  Dw=ceil(abs(normrnd(0.05,0.03))*Nt);
  Pi=randi(c);
  Pl=randi(Nt,1);
  Pw=ceil(abs(normrnd(0.05,0.03))*Nt);
  
  % Randomly select integers to remove from the vector.
  if(Di==1)
      Id=randi(Dw,[1 Nr])+Dl;
  elseif(Di==2)
      %Id=ceil(abs(normrnd(Dl,Dw,[1 Nr])));
      Id=randi(Dw,[1 Nr])+Dl;
  end
  if(max(Id)>Nr)
      Id=Id-(max(Id)-Nr);
  end
  if(min(Id)<1)
      Id=Id+(1-min(Id));
  end
  
  % Remove integers from the vector.
  %for i=1:length(Id)
  %    I(find(I==Id(i),1))=[];
  %end
  I=RepRemove(I,Id);
  
  % Randomly select integers to add to the vector.
  Np=(Nt*Ns)-length(I);
  if(Pi==1)
      Ip=randi(Pw,[1 Np])+Pl;
  elseif(Pi==2)
      Ip=ceil(abs(normrnd(Pl,Pw,[1 Nr])));
  end
  if(max(Ip)>Nr)
      Ip=Ip-(max(Ip)-Nr);
  end
  if(min(Ip)<1)
      Ip=Ip+(1-min(Ip));
  end

  % Add integers to the vector.
  I=[I,Ip];
  
end


% Remove the elements of B from A, while preserving repetitions.
function [A] = RepRemove(A, B)
  while ~isempty(B)
      [~,idx,~] = intersect(A,B);
      A(idx)=[];
      [~,idx,~] = intersect(B,B);
      B(idx)=[];
  end
end

