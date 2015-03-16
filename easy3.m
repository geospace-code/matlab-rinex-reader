% EASY3	  Read RINEX navigation file reformat into Matlab Eph matrix.
%         Open a RINEX observation file analyse the header and identify
%         observation types. The function fepoch_0 finds epoch time 
%         and observed PRNs in an OK epoch (digit 0, RTK observations 
%         will have a 2 in this place). Next we read the observations 
%         and use recpo_ls to get a least-squares estimate for the 
%         (stand alone) receiver position.

%Kai Borre 31-10-2001
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2001/10/31  $

% Read RINEX ephemerides file and convert to
% internal Matlab format
rinexe('site247j.01n','eph.dat');
Eph = get_eph('eph.dat');

% We identify the observation file and open it
ofile1 = 'site247j.01o';
fid1 = fopen(ofile1,'rt');
[Obs_types1, ant_delta1, ifound_types1, eof11] = anheader(ofile1);
NoObs_types1 = size(Obs_types1,2)/2;
Pos = [];

% There are 20 epochs of data in ofile1
for q = 1:20
    [time1, dt1, sats1, eof1] = fepoch_0(fid1);
    NoSv1 = size(sats1,1);
    % We pick the observed P2 pseudoranges
    obs1 = grabdata(fid1, NoSv1, NoObs_types1);
    i = fobs_typ(Obs_types1,'P2');
    pos = recpo_ls(obs1(:,i),sats1,time1,Eph);
    Pos = [Pos pos];
end
me = mean(Pos,2);
display('Mean Position as Computed From 20 Epochs:')
display(['X: ',num2str(me(1,1)),'  Y: ',num2str(me(2,1)),'  Z: ',num2str(me(3,1)) ])
plot((Pos(1:3,:)-Pos(1:3,1)*ones(1,q))','linewidth',2)
title('Positions Over Time','fontsize',16)
legend('X','Y','Z')
xlabel('Epochs [1 s interval]','fontsize',16)
ylabel('Variation in Coordinates, Relative to the First Epoch [m]','fontsize',16)
set(gca,'fontsize',16)
legend