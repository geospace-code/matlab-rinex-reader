%SIGMA_TR Application of Variance propagation.
%         The data are taken from Example 9.2 in
%         Kai Borre (1995): GPS i landmaalingen

%Written by Kai Borre
%Copyright (c) by Kai Borre
%$Revision: 1.0  $  $Date:1998/11/08  $

phi = dms2rad(56, 55, 38.9);
lambda = dms2rad(10, 1, 57.97);
F = [ -sin(lambda) -sin(phi)*cos(lambda) cos(phi)*cos(lambda);
       cos(lambda) -sin(phi)*sin(lambda) cos(phi)*sin(lambda);
                0               cos(phi)             sin(phi)]
Sigma = [  25     -7.97 18.22;
           -7.97   4    -6.36;
           18.22  -6.36 16   ];
F_inv = inv(F);  % because of bad notation in the reference
Sigma_UTM=F_inv*Sigma*F_inv'
std_dev = sqrt(diag(Sigma_UTM))
%%%%%%%%%%%%%%%%%%% end sigma_tr.m  %%%%%%%%%%%%%%%%%%%%