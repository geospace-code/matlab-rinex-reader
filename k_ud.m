function [X,P] = k_ud(X,P,H,Y,R)
%  K_UD   Kalman update, one measurement per call
%	  Observation covariance R

%  Written by Kai Borre and C.C. Goad
%  November 24, 1996

     omc = Y-H'*X;
     HP = H'*P;
     innovation_variance = HP*H+R;
     K = HP'/innovation_variance;
     X = X+K*omc;
     P = P-K*HP;

%%%%%%%% end k_ud.m  %%%%%%%%%%%%%%%%%%
