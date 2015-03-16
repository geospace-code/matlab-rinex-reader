function [omc,bas] = baseline(X_i,obs1,obs2,sats,time,Eph)
% BASELINE Computation of baseline between master and rover 
%          from pseudoranges alone

%Kai Borre 31-10-2001
%Copyright (c) by Kai Borre
%$Revision: 1.1 $  $Date: 2002/11/24  $

m = size(obs1,1);  % number of svs
% identify ephemerides columns in Eph
for t = 1:m
    col_Eph(t) = find_eph(Eph,sats(t),time);
end
% preliminary guess for receiver position 
X_j = X_i(1:3,1);
% Computation of weight matrix
D = [ones(m-1,1) -eye(m-1) -ones(m-1,1) eye(m-1)];
C = inv(D*D');

for iter = 1:8
    % k is the reference satellite. We select the first one
    [tcorr,rhok_j,Xk_ECF] = get_rho(time, obs2(1), Eph(:,col_Eph(1)),X_j);
    [tcorr,rhok_i,Xk_ECF] = get_rho(time, obs1(1), Eph(:,col_Eph(1)),X_i);
    for t = 2:m % t runs over PRNs given in sats; ref.sat. is number 1
        [tcorr,rhol_j,Xl_ECF] = get_rho(time, obs2(t), Eph(:,col_Eph(t)), X_j);
        [tcorr,rhol_i,Xl_ECF] = get_rho(time, obs1(t), Eph(:,col_Eph(t)), X_i); 
        A(t-1,:) = [(Xk_ECF(1)-X_j(1))/rhok_j - (Xl_ECF(1)-X_j(1))/rhol_j,  ...
                    (Xk_ECF(2)-X_j(2))/rhok_j - (Xl_ECF(2)-X_j(2))/rhol_j,  ...
                    (Xk_ECF(3)-X_j(3))/rhok_j - (Xl_ECF(3)-X_j(3))/rhol_j];
        observed = obs1(1)-obs2(1)-obs1(t)+obs2(t);
        calculated = rhok_i-rhok_j-rhol_i+rhol_j;
        omc(t-1,1) = observed - calculated;
    end; % t 
    x = inv(A'*C*A)*A'*C*omc; 
    X_j = X_j+x; 
end % iter    
bas = X_i(1:3,1)-X_j;
%%%%%%%%%%%%%%%%%%%%%  baseline.m  %%%%%%%%%%%%%%%%%%%%%

