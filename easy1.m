%EASY1	  Computation of the essential parameter: 
%         seconds of week, sow

%Kai Borre 27-07-2002
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2002/07/27  $

%Copy of line 29 in the RINEX file site247j.01o
%01  9  4  9 40  0.0000000  0  7G 1G 4G 7G13G20G24G25               
% Compute sow for first epoch in observation file
jd = julday(2001,9,4,9+40/60);
[week,sow] = gps_time(jd)%;
%%%%%%%%%%%%%%%%%%%%% end easy1.m %%%%%%%%%%%%%%%



