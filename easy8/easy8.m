% EASY8k   Test for cycle slip and repair of receiver clock offset 
%          after idea by Kees de Jong.

%Kai Borre 26-12-2002
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2002/12/26  $

v_light = 299792458;    % vacuum speed of light m/s
f0 = 10.23*10^6;
f1 = 154*f0;
f2 = 120*f0;
lambda1 = v_light/f1;
lambda2 = v_light/f2;
alpha = (f1/f2)^2;

ofile1 = 'site24~1.01o';   
fid1 = fopen(ofile1,'rt');
[Obs_types1, ant_delta, ifound_types1, eof11] = anheader(ofile1);
if ((ifound_types1 == 0) | (eof11 == 1))
    error('Basic information is missing in RINEX file')
end;
NoObs_types1 = size(Obs_types1,2)/2;
epochend = 22;

% Next we include the rover observation file and open it
ofile2 = 'SITE247j.01O';
fid2 = fopen(ofile2,'rt');
[Obs_types2, ant_delta2, ifound_types2, eof12] = anheader(ofile2);
NoObs_types2 = size(Obs_types2,2)/2;

j1 = fobs_typ(Obs_types1,'P1');  % pseudorange on L1
k1 = fobs_typ(Obs_types1,'P2');
l1 = fobs_typ(Obs_types1,'L1');  % phase on L1
m1 = fobs_typ(Obs_types1,'L2');

j2 = fobs_typ(Obs_types2,'P1');  % pseudorange on L1
k2 = fobs_typ(Obs_types2,'P2');
l2 = fobs_typ(Obs_types2,'L1');  % phase on L1
m2 = fobs_typ(Obs_types2,'L2');
cols1 = [j1 k1 l1 m1];
cols2 = [j2 k2 l2 m2];

% Preparations for filter
x_1 = zeros(4,1);      % state vector: [B1 B2 B3 Idot]
delta_t = 1;         % epoch interval in seconds
F = eye(4);
F(1:3,4) = delta_t;
A = diag([alpha-1, -2, -alpha-1]); 
A = [A zeros(3,1)];
T = [-ones(3,1) eye(3)];
Sigma_b = 2*diag([0.3^2 0.3^2 0.003^2 0.003^2]);
Sigma_e = T*Sigma_b*T';
Sigma_epsilon = diag([.1^2 .1^2 .1^2 .1^2]);
P_1 = 10^2*eye(4);
X = [];

fighdl = figure; 
set(gcf,'UserData',zeros(4,1)*inf);
pl_handle = plot(1,zeros(4,1)*inf,'.','Erasemode','none','Markersize',5);
axis([0 epochend+1 -7 18]);

for epoch = 1:epochend
    [time1, dt1, sats1, eof1] = fepoch_0(fid1);
    [time2, dt2, sats2, eof2] = fepoch_0(fid2);
    if time1 ~= time2
        disp('Epochs not corresponding')
        break
    end;
    NoSv1 = size(sats1,1);
    NoSv2 = size(sats2,1);
    % We pick the observations
    obs1 = grabdata(fid1, NoSv1, NoObs_types1);
    obs2 = grabdata(fid2, NoSv2, NoObs_types2);
    % Rearranging the obs2 row sequence to correspond to the obs1 row sequence
    obs = obs2;
    for s = 1:NoSv1
        ind = find(sats1(s) == sats2(:));
        obs2(s,1) = obs(ind,1);
    end
    Obs1 = obs1(:,cols1); % selecting and ordering the columns used
    Obs2 = obs2(:,cols2);    
    b = [Obs2(:,2)-Obs2(:,1)-(Obs1(:,2)-Obs1(:,1));                  % P2-P1
        Obs2(:,3)*lambda1-Obs2(:,1)-(Obs1(:,3)*lambda1-Obs1(:,1));   % Phi1-P1
        Obs2(:,4)*lambda2-Obs2(:,1)-(Obs1(:,4)*lambda2-Obs1(:,1))];  % Phi2-P1
    
    % for s = 1:NoSv1
    % We start with the first satellite
    s = 1; %%%%%%%%%%%%%%%%%
    b_1 = [b(s);b(s+NoSv1);b(s+2*NoSv1)];
    % Kalman filter, first filtering
    PAt = P_1*A';
    Ivar = A*PAt+Sigma_e;
    K_1 = PAt*inv(Ivar);
    x_1 = x_1+K_1*(b_1-A*x_1);
    P_1 = P_1-K_1*A*P_1;
    % next prediction
    x_1 = F*x_1;
    X = [X x_1];
    P_1 = F*P_1*F'+Sigma_epsilon;
    set(pl_handle,'xdata',epoch*ones(4,1),'ydata',X(:,epoch)');
    drawnow   
    % end % s   
end
fclose all;

ylabel('Changes in {\itB}_1, {\itB}_2, {\itB}_3, and {\itdI/dt} [m]')
xlabel('Epochs, [1 s interval]')

fighdl2 = figure;
plot(1:epoch,X','linewidth',2) % ,'-'
title(['Check of Cycle Slips for PRN  ' num2str(sats1(s))],'fontsize',16)
ylabel('Variations in {\itB}_1, {\itB}_2, {\itB}_3, and {\itdI/dt}  [m]','fontsize',16)
xlabel('Epochs [1 s interval]','fontsize',16)
legend('{\itB}_1','{\itB}_2','{\itB}_3','{\itdI/dt}') 
set(gca,'fontsize',16)
legend

print -deps easy8

break
% Repair of clock reset of 1ms ~ 299 km; affects only pseudoranges
i1 = find(abs(DP(1,:)) > 280000);

for j = i1
    if DP(:,j) < 0
        DP(:,j) = DP(:,j)+299792.458; 
    else 
        DP(:,j) = DP(:,j)-299792.458; 
    end
end

figure(1);
%colorordermatrix = [.97 .4  .25;
%    .15 .25 .09; 
%    .29 .25 .09;
%    .55 .7  .51;
%    .23 .73 1;
%    .31 .57 .35;
%    .72 .68 .35;
%    .99 .93 .96;
%    .51 .02 .25;
%    .08 .11 .33;
%    .73 .23 .56;
%    .97 .33 .19];
%axes('Colororder',colorordermatrix,'NextPlot','add');
plot((deltaP-deltaPhi)','linewidth',2)
legend(eval('num2str(sv)'),2)
title('Check of Cycle Slips')
ylabel('Misclosure [m]')
xlabel('Epochs  [Interval 1 s]')
legend

%%%%%%%%% end easy8k.m %%%%%%%%%

