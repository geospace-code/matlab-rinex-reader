function bessel_1(phi1d,phi1m,phi1s,l1d,l1m,...
                                    l1s,A1d,A1m,A1s,s12,a,finv)
%BESSEL_1 Solution of the direct geodetic problem according to
%         the Bessel-Helmert method as described in Zhang Xue-Lian.
%         Given a point with coordinates (phi1, l1) and 
%         a gedesic with azimuth A1 and length s12 from here. The 
%         given reference ellipsoid has semi-major axis a and 
%         inverse flattening finv. The coordinates and the azimuth
%         have the format degree, minute, and second
%         with decimals. 
%
%         Zhan Xue-Lian (1985): The nested coefficient method
%            for accurate solutions of direct and inverse geodetic
%            problems with any length. In proceedings of the 7th
%            International Symposium on Geodetic Computations, 
%            Cracow, June 18--21, 1985. Pages 747--763. Institute
%            of Geodesy and Cartography. Wasaw, Poland, ul. Jasna 2/4.
%
%         A good decription of the general background of the problems
%         is given in
%
%  	    Bodem\"uller, H.(1954): Die geod\"atischen Linien des 
%	         Rotationsellipsoides und die L\"osung der geod\"atischen
%	         Hauptaufgaben f\"ur gro\ss{}e Strecken unter 
%            besonderer Ber\"ucksichtigung der Bessel-Helmertschen
%            L\"osungsmethode. Deutsche Geod\"atische Kommission,
%            Reihe B, Nr. 13.

%Kai Borre, January 24, 1999
%Copyright (c) by Kai Borre
%$Revision 1.1 $  $Date:2002/05/02  $

dtr = pi/180;        % degrees to radians
if nargin == 0
phi1  = 50*dtr;
l1 =    10*dtr;
A1 =   140*dtr;
s12 = 15000000;   % m
a = 6378388;
finv = 297;
else
   phi1 = dms2rad(phi1d,phi1m,phi1s);
   l1 = dms2rad(l1d,l1m,l1s);
   A1 = dms2rad(A1d,A1m,A1s);
end

f = 1/finv; 
ex2 = (2-f)*f/(1-f)^2;	        % second eccentricity squared
tanu1 = (1-f)*tan(phi1);        % (1)
sigma1 = atan2(tanu1,cos(A1));  % (2)
u1 = atan(tanu1);
cosun = cos(u1)*sin(A1);        % (3)
sinun2 = 1-cosun^2;
t = ex2*sinun2/4;               % (4)
K1 = 1+t*(1-t*(3-t*(5-11*t))/4);
K2 = t*(1-t*(2-t*(37-94*t)/8));
v = f*sinun2/4;                 % (5)
K3 = v*(1+f+f^2-v*(3+7*f-13*v));
deltasigma_old = 0;
deltasigma_new = 1;

while  abs(deltasigma_old-deltasigma_new) > 1.e-12
   deltasigma_old = deltasigma_new;
   sigma = s12/(K1*(1-f)*a)+deltasigma_old; % (6)
   sigmam = 2*sigma1+sigma;
   deltasigma_new = K2*sin(sigma)*(cos(sigmam)+...
                      K2*(cos(sigma)*cos(2*sigmam)+...
                        K2*(1+2*cos(2*sigma))*cos(3*sigmam)/6)/4); %(7)
end

tanu2 = (sin(u1)*cos(sigma)+cos(u1)*sin(sigma)*cos(A1))/...
                           sqrt(1-sinun2*(sin(sigma1+sigma))^2); %(8)
phi2 = atan(tanu2/(1-f));
disp('Phi2');
rad2dms(phi2);
deltaomega = (1-K3)*f*cosun*(sigma+...
               K3*sin(sigma)*(cos(sigmam)+...
               K3*cos(sigma)*cos(2*sigmam))); % (9)
omega = atan2(sin(sigma)*sin(A1),...
               cos(u1)*cos(sigma)-sin(u1)*sin(sigma)*cos(A1)); % (10)
l2 = l1+omega-deltaomega;
disp('lambda2')
rad2dms(l2);
A2 = atan2(cos(u1)*sin(A1),...
               cos(u1)*cos(sigma)*cos(A1)-sin(u1)*sin(sigma)); % (11)             
disp('A2');
rad2dms(A2);              

%----------------------------------------------

function result = dms2rad(deg,min,sec);
% Conversion of degrees, minutes, and seconds to radians

neg_arg = 'FALSE';
if deg < 0
   neg_arg = 'TRUE ';
   deg = -deg;
end
arg = deg+min/60+sec/3600;
result = arg*pi/180;
if neg_arg == 'TRUE ';
   result = -result;
end

%------------------------------------------

function result = rad2dms(arg)
%RAD2DMS Conversion of radians to degrees, minutes, and seconds%

neg_arg = 'FALSE';
if arg < 0
   neg_arg = 'TRUE ';
   arg = -arg;
end

arg = arg*180/pi;
result = zeros(1,3);
result(1) = fix(arg);
if result(1) == 0
   result(2) = fix(arg*60);
else
   result(2) = fix(rem(arg,result(1))*60);
end
result(3) = (arg-result(1)-result(2)/60)*3600;
if neg_arg == 'TRUE '
   result(1) = -result(1);
end
fprintf('   %3.0f %2.0f %8.6f\n',result(1),result(2),result(3))

%%%%%%%%%%%%%%%%% end bessel_1.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%

