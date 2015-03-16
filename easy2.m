%EASY2	 Convert observation time into sow.
%        We read the corresponding RINEX navigation file 
%        and reformat the data into the Matlab matrix Eph.
%        For given SV we find the corresponding column in Eph 
%        and call the basic satpos function

%Kai Borre 27-07-2002
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2002/07/27  $

% Compute sow for first epoch in observation file site247j.01o
% 01  9  4  9 40  0.0000000  0  7G 1G 4G 7G13G20G24G25               
jd = julday(2001,9,4,9+40/60);
[week,sow] = gps_time(jd);

% Read RINEX ephemerides file and convert to
% internal Matlab format
rinexe('SITE247J.01N','eph.dat');
Eph = get_eph('eph.dat');

% We identify the observed satellites in line 29 of RINEX file site247j.01o
% 01  9  4  9 40  0.0000000  0  7G 1G 4G 7G13G20G24G25               
svs = [1 4 7 13 20 24 25];
for t = 1:length(svs)
    col_Eph(t) = find_eph(Eph,svs(t),sow);
    sat(1:3,t) = satpos(sow,Eph(:,col_Eph(t)));
end

sat     % position of svs in ECEF system
%%%%%%%%%%%%%%%%%%%%% end easy2.m %%%%%%%%%%%%%%%



