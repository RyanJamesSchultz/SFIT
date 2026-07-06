function D=PreProcData(Cases,Si)
  % Function that will pre-process the catalogue and injection information.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Get the data table and define an output strucutre.
  table_CaseDataTable;
  D=struct('Case',[],'Stages',[],'Lat',[],'Lon',[],'Dep',[],'T',[],'M',[],'Mc',[],'Mk',[],'t',[],'v',[],'V',[],'T1',[],'T2',[],'T3',[]);
  
  % Loop over all the user-input cases.
  for i=1:length(Cases)
      
      % Get the table entry for this case.
      Case=Cases{i};
      D(i).Case=Case;
      D(i).Stages=Si;
      j=find(strcmpi(Case,{T.Name}));
      load(T(j).Filename);
      D(i).Mc=T(j).Mc;
      D(i).Mk=T(j).Mk;
      D(i).dMb=T(j).dMb;
      D(i).dMd=T(j).dMd;
      
      % Catalogue information.
      D(i).Lat=S.cat.lat;
      D(i).Lon=S.cat.lon;
      D(i).Dep=S.cat.dep;
      D(i).T=S.cat.time;
      D(i).M=S.cat.mag;

      % Dither the times.
      dtd=T(j).dtd;
      f=24*3600;
      D(i).T=D(i).T+dtd*rand(size(D(i).T))/f-(dtd/2)/f;

      % Dither the magnitudes.
      dMd=T(j).dMd;
      D(i).M=D(i).M+(dMd)*rand(size(D(i).M))-dMd/2;
      
      % Chronologically sort the catalogue.
      [~,I]=sort(D(i).T);
      D(i).Lat=D(i).Lat(I);
      D(i).Lon=D(i).Lon(I);
      D(i).Dep=D(i).Dep(I);
      D(i).T=D(i).T(I);
      D(i).M=D(i).M(I);
      
      % Well trajectories.
      D(i).Wlat=S.traj.lat;
      D(i).Wlon=S.traj.lon;
      D(i).Wtvd=S.traj.dep;
      D(i).Wmd=S.traj.md;
      
      % Stage MD locations.
      D(i).MDs=T(j).md(Si);
      
      % Injection information.
      D(i).t=S.inj.time; % days.
      D(i).v=S.inj.rate; % m³/min.
      
      % Convert to new time format (if necessary).
      if(~isdatetime(D(i).T))
          D(i).T=datetime(D(i).T,'ConvertFrom','datenum');
          D(i).t=datetime(D(i).t,'ConvertFrom','datenum');
      end
      
      % Get the start and end time of the stages.
      D(i).T1=min(T(j).Ts(Si));
      D(i).T2=max(T(j).Te(Si));
      D(i).Tf=max(T(j).Tf);
      
      % Get the end time of the filter window.
      k=max(Si);
      if(k==length(T(j).Ts))
          D(i).T3=max(D(i).Tf);
      else
          D(i).T3=min([T(j).Ts(k+1) max(D(i).Tf)]);
      end
      
      % Temporally filter the catalogue over the appropriate timeframe.
      I=(D(i).T>=D(i).T1)&(D(i).T<=D(i).T3);
      D(i).Lat=D(i).Lat(I);
      D(i).Lon=D(i).Lon(I);
      D(i).Dep=D(i).Dep(I);
      D(i).T=D(i).T(I);
      D(i).M=D(i).M(I);
      
      % Temporally filter the injection information over the appropriate timeframe.
      I=(D(i).t>=D(i).T1)&(D(i).t<=D(i).T3);
      D(i).t=D(i).t(I);
      D(i).v=D(i).v(I);
      
      % Compute cumulative volume (m³).
      f=60*24; % minutes/day
      dt=days(D(i).t(2)-D(i).t(1));
      D(i).V=cumsum(D(i).v)*dt*f;
      D(i).v=D(i).v*f; % m³/day.
      
  end
  
end


