% EASY10 plots ionospheric delay from phase observations.
%        The present code does not handle
%   	 1. cycle slips, and
%	     2. outliers.

%Kai Borre 27-07-2002
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2002/07/27  $

% Initial computations of constants
v_light = 299792458;	     % vacuum speed of light m/s
f1 = 154*10.23E6;		     % L1 frequency Hz
f2 = 120*10.23E6;			 % L2 frequency Hz
lambda1 = v_light/f1;	     % wavelength on L1:  .19029367  m
lambda2 = v_light/f2;	     % wavelength on L2:  .244210213 m

% We identify the observation file and open it
ofile = 'kofi1.01o'; 
fid = fopen(ofile,'rt');
[Obs_types, ant_delta, ifound_types, eof1] = anheader(ofile);
NoObs_types = size(Obs_types,2)/2;
obsstr(1,1:2) = 'P1'; % P1
obsstr(2,1:2) = 'P2'; % P2
obsstr(3,1:2) = 'L1'; % Phi1
obsstr(4,1:2) = 'L2'; % Phi2
match = zeros(1,4);
for t = 1:4
    for ii = 1:NoObs_types
        mat = strmatch(obsstr(t,1:2),Obs_types(1,2*ii-1:2*ii),'exact');
        if isempty(mat) == 0, match(1,t) = ii; end
    end
end
Oc = match; 
qend = 300;
phase_diff = zeros(32,qend)*inf;
for q = 1:qend
    [time, dt, sats, eof] = fepoch_0(fid);
    NoSv = size(sats,1);
    obsm = grabdata(fid, NoSv, NoObs_types);
    obs = obsm(:,Oc); % P1 P2 Phi1 Phi2
    for t = 1:NoSv
        Phi1 = obs(t,3)*lambda1; 
        Phi2 = obs(t,4)*lambda2; 
        sat = sats(t);
        phase_diff(sat,q) = Phi1-Phi2;
    end
end; 
pd = phase_diff-phase_diff(:,1)*ones(1,qend);
prn = find(~isnan(pd(:,1)));
% We must delete rows with NaN's in order to cycle correctly through the colours
pd(isnan(pd(:,1:qend))) = []; 
pd = reshape(pd,length(prn),qend);
pd = pd/(1-(f1/f2)^2);

figure(1);
%colorordermatrix = [.97 .4  .25;
%                    .15 .25 .09; 
%                    .29 .25 .09;
%                    .55 .7  .51;
%                    .23 .73 1;
%                    .31 .57 .35;
%                    .72 .68 .35;
%                    .99 .93 .96;
%                    .51 .02 .25;
%                    .08 .11 .33;
%                    .73 .23 .56;
%                    .97 .33 .19];
%axes('Colororder',colorordermatrix,'NextPlot','add');
plot(1000*pd','linewidth',2)
legend(eval('num2str(prn)'),2)
title('Ionospheric Delay From {\itL}_1 and {\itL}_2 Phase','fontsize',16)
ylabel('Ionospheric Delay  [mm]','fontsize',16)
xlabel('Epochs  [1 s interval]','fontsize',16)
set(gca,'fontsize',16)
legend
print -deps easy101

%figure(2);
%plot(1:300,fft(1000*pd'))
for qq = 1:7
    autocorr(1000*pd(qq,:)');
end
%print -deps easy102

%%%%%%%%%%%%%%%%%%%%%% end easy10.m  %%%%%%%%%%%%%%%%%%%







