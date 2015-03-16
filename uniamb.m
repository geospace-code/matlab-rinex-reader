%Script to do examples in 
% L.T. Liu, H.T. Hsu, Y.Z. Zhu, J.K. Ou (1999):
% A new approach to GPS ambiguity decorrelation
% Journal of Geodesy, 73: 478---490

% Copyright (c) by Kai Borre
% 17-Nov-1999, Helsinki

% Example 1

Qa = [100 19; 19 4];
[l,u,p] = lu(Qa);
lr = round(l);
Qh = lr*Qa*lr';


%Example 2

Qa = [100 19; 19 4];
L = eye(size(Qa,1));
L([1 2],1:2) = L([2 1],1:2);
Ql = L*Qa*L';
[H,U] = slu(Ql);
hr = round(H);
Qu = hr*Ql*hr';


%Example 3

Qa = [80.00  9.75 27.80;
   9.75  1.25  3.45;
   27.80  3.45 10.00];
n = size(Qa,1);
dia = diag(Qa);
D = Qa;
L_acc = [];
K_acc = [];

%Step 1: Permutation of covariance matrix
for i = 1:n
   [k0,j] = sort(dia);
   L = eye(n);
   if i < j(1)
      L(:,[i j(1)]) = L(:,[j(1) i]);
   end
   Ql = L*D*L';
   L_acc = [L_acc L];
   K = eye(n);
   K(i+1:n,i) = -Ql(i+1:n,i)/Ql(i,i);
   K_acc = [K_acc K];
   D = K*Ql*K';
   dia = diag(D(i:n,i:n));
end
D

%Step 2: Float decorrelation

H = eye(n);
for i = 1:n
   H = H*L_acc(:,3*(i-1)+1:3*i)*K_acc(:,3*(i-1)+1:3*i);
end
H

%Step 3: United Ambiguity decorrelation


%%%%%%%%%%%%%%%%%% end uniamb.m  %%%%%%%%%%%%%%%%%%