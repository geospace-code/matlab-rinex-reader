function [x,P] = k_updatm(x,P,A,b,R,Q)
%K_UPDATM Kalman update for matrix
%	       Allows for system covariance Q
%	       Allows for observation covariance R

%Kai Borre, February 21, 1999

omc = b-A*x;
P = P + Q;
PAt = P*A';
Ivar = A*PAt+R;
K = PAt*inv(Ivar);
x = x+K*omc;
P = P-K*A*P;
%%%%%%%% end k_updatm.m  %%%%%%%%%%%%%%%%%%
