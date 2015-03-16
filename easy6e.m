%EASY6E	We compute a baseline from C/A code and phase observations.
%	    The code does not handle
%	       1. cycle slips and
%	       2. outliers.
%	    In contrast to EASY5, we now use an extended Kalman filter for the
%	    estimation. The present code is no real RTK code as all computational
%	    steps do not happen on an epoch-by-epoch basis
%
%       The filter uses e(ponentially) age-weighting of old data.

%       Reference:
%       Kailath, Thomas & Ali H. Sayed & Babak Hassibi (2000): Linear Estimation. 
%       Prentice Hall. Pages 68--69

%Kai Borre 27-07-2002
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2002/07/27  $

% Initial computations of constants
v_light = 299792458;	     % vacuum speed of light m/s
f1 = 154*10.23E6;		     % L1 frequency Hz
f2 = 120*10.23E6;			 % L2 frequency Hz
lambda1 = v_light/f1;	     % wavelength on L1:  .19029367  m
lambda2 = v_light/f2;	     % wavelength on L2:  .244210213 m

% Read RINEX ephemerides file and convert to internal Matlab format
rinexe('SITE247J.01N','eph.dat');
Eph = get_eph('eph.dat');

% We identify the master observation file and open it
ofile1 = 'SITE247J.01O';
fid1 = fopen(ofile1,'rt');
[Obs_types1, ant_delta1, ifound_types1, eof11] = anheader(ofile1);
NoObs_types1 = size(Obs_types1,2)/2;

% We start by estimating the master position
[time1, dt1, sats1, eof1] = fepoch_0(fid1);
NoSv1 = size(sats1, 1);
m = NoSv1;
obs1raw = grabdata(fid1, NoSv1, NoObs_types1);
i = fobs_typ(Obs_types1,'C1'); % We use C/A pseudoranges
[X_i, el] = recpo_ls(obs1raw(:,i), sats1, time1, Eph);
[phi_i,lambda_i,h_i] = togeod(6378137,298.257223563,X_i(1),X_i(2),X_i(3));
% We close all files to ensure that the next reading starts
% at the top of the observation files
fclose all;

% Finding columns in Eph for each SV
for t = 1:m
    col_Eph(t) = find_eph(Eph,sats1(t),time1);
end

% Computation of elevation angle to all SVs.
all_sats1 = sats1;
% Delete Sv with elevation smaller than 10 degrees
sats1(el<10) = [];
del_sat = setdiff(all_sats1,sats1);

no_del_sat = [];
for t = 1:length(del_sat)
    no_dels = find(del_sat(t) == all_sats1);
    no_del_sat = [no_del_sat; no_dels];
end
No_del_sat = length(no_del_sat);

% Selecting reference SV. We take the SV with largest elevation
[y,ind] = max(el);
refsv = sats1(ind);
ofile1 = 'SITE247J.01O';
fid1 = fopen(ofile1,'rt');
ofile2 = 'SITE24~1.01O';
fid2 = fopen(ofile2,'rt');

% We start reading both observation files
[Obs_types1, ant_delta1, ifound_types1, eof11] = anheader(ofile1);
NoObs_types1 = size(Obs_types1,2)/2;
obsstr = ['P1';'P2';'L1';'L2'];
Oc = [];
for t = 1:4
    oc = strmatch(obsstr(t,:),strvcat('C1','P1','P2','L1','L2'),'exact');
    Oc = [Oc oc];
end
[Obs_types2, ant_delta2, ifound_types2, eof12] = anheader(ofile2);
NoObs_types2 = size(Obs_types2,2)/2;

% Computation of covariance matrix Sigma for double differenced observations
m1 = m-No_del_sat; % original number of SVs - deleted SVs due to low elevations
D = [ones(m1,1) -eye(m1) -ones(m1,1) eye(m1)];
Sigma = D*D';
X = zeros(3+2*m1,1);	      % coord.diff., N1, N2
N = zeros(3+2*m1,3+2*m1);     % initialization of normals
rs = zeros(3+2*m1,1);	      % initialization of right side
X_a = [];
X_j = X_i(1:3,1);
refrow = find(refsv == sats1);

% We process three epochs for estimating ambiguities; the present data evidently
% need three or more epochs for getting reliable estimates of the float ambiguities
for q = 1:6
    X_j = X_i(1:3,1)+X(1:3,1);
    [time1, dt1, sats1, eof1] = fepoch_0(fid1);
    [time2, dt2, sats2, eof2] = fepoch_0(fid2);
    if time1 ~= time2
	disp('Epochs do not correspond in time')
	break
    end;
    time = time1;
    NoSv1 = size(sats1,1);
    NoSv2 = size(sats2,1);
    obsm = grabdata(fid1, NoSv1, NoObs_types1);
    obsr = grabdata(fid2, NoSv2, NoObs_types2);
    obs1 = obsm(:,Oc); % P1 P2 Phi1 Phi2
    % Reordering of rows in obsr to correspond to obsm
    for s = 1:NoSv1
	Ind = find(sats1(s) == sats2(:));
	obs2(s,:) = obsr(Ind,Oc);
    end
    % Computing rho for refsv
    [tcorr,rhok_j,Xk_ECF] = get_rho(time, obs2(refrow,1), Eph(:,col_Eph(refrow)), X_j);
    [tcorr,rhok_i,Xk_ECF] = get_rho(time, obs1(refrow,1), Eph(:,col_Eph(refrow)), X_i);
    tt = 0;
    A1 = [];
    t0 = 1:NoSv1;
    t1 = setdiff(t0,no_del_sat); % we delete the low satellites
    for t = t1
	tt = tt+1;
	[tcorr,rhol_j,Xl_ECF] = get_rho(time,obs2(t,1), Eph(:,col_Eph(t)), X_j);
	[tcorr,rhol_i,Xl_ECF] = get_rho(time,obs1(t,1), Eph(:,col_Eph(t)), X_i);
	A0 = [(Xk_ECF(1)-X_j(1))/rhok_j - (Xl_ECF(1)-X_j(1))/rhol_j  ...
		(Xk_ECF(2)-X_j(2))/rhok_j - (Xl_ECF(2)-X_j(2))/rhol_j ...
		(Xk_ECF(3)-X_j(3))/rhok_j - (Xl_ECF(3)-X_j(3))/rhol_j];
	A1 = [A1; A0];
	Phi1 = (obs1(refrow,3)-obs1(t,3)-obs2(refrow,3)+obs2(t,3))*lambda1;
	Phi2 = (obs1(refrow,4)-obs1(t,4)-obs2(refrow,4)+obs2(t,4))*lambda2;
	b(tt,:) = Phi1-lambda1*X(3+tt,1);
	b(m1+tt,:) = Phi2-lambda2*X(3+m1+tt,1);
	bk(tt,:) =  rhok_i-rhok_j-rhol_i+rhol_j;
	bk(m1+tt,:) =  rhok_i-rhok_j-rhol_i+rhol_j;
    end;
    A_modi = eye(m1);		 % modified coefficient matrix
    col = find(refsv == sats1);  % find column for reference PRN
    A_modi(:,col) = -ones(m1,1);
    A_aug = [A1 lambda1*A_modi 0*eye(m1); A1 0*eye(m1) lambda2*A_modi];
    N = N+A_aug'*kron(eye(2),Sigma)*A_aug;
    rs = rs+A_aug'*kron(eye(2),Sigma)*(b-bk);
end %q
PP = pinv(N);
% X contains the three preliminary baseline components and the float ambiguities
X = PP*rs;

% Estimation of ambiguities by means of the Lambda method
[a,sqnorm,Sigma_afixed,Z] = lambda(X(4:4+2*m1-1,1),PP(4:4+2*m1-1,4:4+2*m1-1));
% Correcting baseline vector as consequence of changing float ambiguities to fixed ones
X(1:3,1) = X(1:3,1)-PP(1:3,4:4+2*m1-1)*inv(PP(4:4+2*m1-1,4:4+2*m1-1))*...
			     (X(4:4+2*m1-1,1)-a(:,1)); %select first set of candidates
X(4:4+2*m1-1,1) = a(:,1);
fprintf('\n N1 for PRN %3.0f: %3.0f',[sats1(t1)'; a(1:m1,1)'])
fprintf('\n')
fprintf('\n N2 for PRN %3.0f: %3.0f',[sats1(t1)'; a(m1+1:2*m1,1)'])

% We close and reopen all files in order to start reading at a known position
fclose all;
ofile1 = 'SITE247J.01O';
fid1 = fopen(ofile1,'rt');
ofile2 = 'SITE24~1.01O';
fid2 = fopen(ofile2,'rt');

% Setting covariances for the Kalman filter; the state vector contains (x,y,z)
P = eye(3);	                    		 % covariances of state vector
Q = 0.05^2*eye(3);			             % covariances of system
R = 0.005^2*kron(eye(2),inv(Sigma));	 % covariances of observations
% In ofile2 we substitute empty observations with NaN's to obtain 22 valid epochs
qend = 22;
% Preliminary estimate of baseline components
x = X(1:3,1);
x_acc = [];
delta_x = zeros(3,1);

for q = 1:qend
    X_j = X_i(1:3,1)+x;
    [phi_j,lambda_j,h_j] = togeod(6378137,298.257223563,X_j(1),X_j(2),X_j(3));
    [time1, dt1, sats1, eof1] = fepoch_0(fid1);
    [time2, dt2, sats2, eof2] = fepoch_0(fid2);
    if time1 ~= time2
	disp('Epochs do not correspond in time')
	break
    end;
    time = time1;
    NoSv1 = size(sats1,1);
    NoSv2 = size(sats2,1);
    obsm = grabdata(fid1, NoSv1, NoObs_types1);
    obsr = grabdata(fid2, NoSv2, NoObs_types2);
    obs1 = obsm(:,Oc); % P1 P2 Phi1 Phi2
    % Reordering of rows in obsr to correspond to obsm
    for s = 1:m
	Ind = find(sats1(s) == sats2(:));
	obs2(s,:) = obsr(Ind,Oc);
    end
    % Computing rho for refsv
    [tcorr,rhok_j,Xk_ECF] = get_rho(time, obs2(1,1), Eph(:,col_Eph(1)), X_j);
    [tcorr,rhok_i,Xk_ECF] = get_rho(time, obs1(1,1), Eph(:,col_Eph(1)), X_i);
    tt = 0;
    A = zeros(2*m1,3);
    for t = t1
	tt = tt+1;
	[tcorr,rhol_j,Xl_ECF] = get_rho(time,obs2(t,1), Eph(:,col_Eph(t)), X_j);
	[tcorr,rhol_i,Xl_ECF] = get_rho(time,obs1(t,1), Eph(:,col_Eph(t)), X_i);
	A0 = [(Xk_ECF(1)-X_j(1))/rhok_j - (Xl_ECF(1)-X_j(1))/rhol_j  ...
		(Xk_ECF(2)-X_j(2))/rhok_j - (Xl_ECF(2)-X_j(2))/rhol_j ...
		(Xk_ECF(3)-X_j(3))/rhok_j - (Xl_ECF(3)-X_j(3))/rhol_j];
	A(tt,:) =  A0;
	A(m1+tt,:) = A0;
	% Tropospheric correction using standard meteorological parameters
	%[az,el_ki,d] = topocent(X_i(1:3),Xk_ECF-X_i(1:3));
	%[az,el_li,d] = topocent(X_i(1:3),Xl_ECF-X_i(1:3));
	%[az,el_kj,d] = topocent(X_j(1:3),Xk_ECF-X_j(1:3));
	%[az,el_lj,d] = topocent(X_j(1:3),Xl_ECF-X_j(1:3));
	%el_ki,    el_li,    el_kj,    el_lj
	%t_corr = tropo(sin(el_lj*pi/180),...
	%    h_j*1.e-3,1013,293,50,0,0,0)...
	%    -tropo(sin(el_li*pi/180),....
	%    h_i*1.e-3,1013,293,50,0,0,0)...
	%    -tropo(sin(el_kj*pi/180),...
	%    h_j*1.e-3,1013,293,50,0,0,0)...
	%    +tropo(sin(el_ki*pi/180),...
	%    h_i*1.e-3,1013,293,50,0,0,0);
	Phi1 = (obs1(refrow,3)-obs1(t,3)-obs2(refrow,3)+obs2(t,3))*lambda1; %-t_corr;
	Phi2 = (obs1(refrow,4)-obs1(t,4)-obs2(refrow,4)+obs2(t,4))*lambda2; %-t_corr;
	b(tt,:) = Phi1-lambda1*a(tt,1);
	b(m1+tt,:) = Phi2-lambda2*a(m1+tt,1);
	bk(tt,:) =  rhok_i-rhok_j-rhol_i+rhol_j;
	bk(m1+tt,:) =  rhok_i-rhok_j-rhol_i+rhol_j;
    end; % t
    
    % Age weighting of old data, see Fagin.
    tau = .1; % The smaller tau is, the faster old observations are forgotten
    age_weight = exp(1/tau);  
    
    %Extended Kalman filter, see pages 509--510 in Strang & Borre (1997): Linear
    % Algebra, Geodesy, and GPS, Wellesley-Cambridge Press
    P = P+Q;
    K = P*A'*inv(A*P*A'/age_weight+ R)/age_weight;
    x = x+K*(b-bk);
    P = (eye(3)-K*A)*P/age_weight;
    fprintf('\nx: %8.3f m,  y: %8.3f m,  z: %8.3f m\n',x(1),x(2),x(3))
    x_acc = [x_acc x];
end %q

% Transformation of geocentric baseline coordinates into topocentric coordinates
for i = 1:qend
    [e(i),n(i),u(i)] = xyz2enu(phi_j,lambda_j,x_acc(1,i),x_acc(2,i),x_acc(3,i));
end
fprintf('\n\nBaseline Components\n')
fprintf('\nX: %8.3f m,  Y: %8.3f m,  Z: %8.3f m\n', ...
				  x_acc(1,qend),x_acc(2,qend),x_acc(3,qend))
fprintf('\nE: %8.3f m,  N: %8.3f m,  U: %8.3f m\n',mean(e),mean(n),mean(u))

figure;
plot(1:qend,[(e-e(1))' (n-n(1))' (u-u(1))']*1000,'linewidth',2)
title('Estimates of Baseline Using Age-Weighting','fontsize',16)
ylabel('State Vector, Changes Relative to Initial Epoch [mm]','fontsize',16)
xlabel('Epochs [1 s interval]','fontsize',16)
legend('East','North','Up')
set(gca,'fontsize',16)
legend

print -deps easy6e
%%%%%%%%%%%%%%%%%%%%%% end easy6e.m  %%%%%%%%%%%%%%%%%%%







