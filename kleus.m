function x2 = kleus(obs)
% KLEUS    Explicit method for computing preliminary receiver
%          coordinates from pseudoranges and ECEF coordinates 
%          of four or more satellites

%  Reference: Kleusberg, A. (1994): Die direkte L\"osung des
%		            	    r\"aumlichen Hyperbelschnitts. Zeitschrift
%			                f\"ur Vermessungswesen, pp. 188--192

% Written by Kai Borre
% May 20, 1997

differ = obs(:,4)-obs(1,4);
d = differ(2:4,1);

b = [  norm(obs(2,1:3)-obs(1,1:3));
       norm(obs(3,1:3)-obs(1,1:3));
       norm(obs(4,1:3)-obs(1,1:3))];
e1 = (obs(2,1:3)-obs(1,1:3))'/b(1);
e2 = (obs(3,1:3)-obs(1,1:3))'/b(2);
e3 = (obs(4,1:3)-obs(1,1:3))'/b(3);
F1 = b(1)/(b(1)^2-d(1)^2)*e1 - b(2)/(b(2)^2-d(2)^2)*e2;
F2 = b(2)/(b(2)^2-d(2)^2)*e2 - b(3)/(b(3)^2-d(3)^2)*e3;
f1 = F1/norm(F1);
f2 = F2/norm(F2);

A = [f1'; f2'];
u = [ (d(2)/(b(2)^2-d(2)^2)-d(1)/(b(1)^2-d(1)^2))/norm(F1);
      (d(3)/(b(3)^2-d(3)^2)-d(2)/(b(2)^2-d(2)^2))/norm(F2)];
g = cross(f1,f2);
h = u(2)*f1-u(1)*f2;
E1 = (cross(g,h) + g*sqrt((norm(g))^2 - (norm(h))^2))/(norm(g))^2;
E2 = (cross(g,h) - g*sqrt((norm(g))^2 - (norm(h))^2))/(norm(g))^2;
s1 = (b(1)^2-d(1)^2)/(2*(d(1)+b(1)*sum(E1.*e1)));
s2 = (b(1)^2-d(1)^2)/(2*(d(1)+b(1)*sum(E2.*e1)));

x2 = obs(1,1:3)'+s1*E1;  % s2
x4 = obs(1,1:3)'+s1*E2;  % s2
%%%%%%%%%%% end kleus.m  %%%%%%%%%%%%%%%%%%%%%%
