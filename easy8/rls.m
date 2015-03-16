function rls(A,b,Sigma)
%RLS      Recursive Least Squares
%         A is the coefficient matrix, b the observations and 
%         Sigma a vector containing the diagonal entries of
%         the covariance matrix for the problem.
%         We include one additional observation for increasing i by 1

%Copyright (c) by Kai Borre
%$Revision: 1.0 $  %Date:1999/10/24  $

if nargin == 0
   A = [1 1;1 2;-1 1];
   b = [2;1;0];
   Sigma = diag([1,.5,1]);
end

% Initial weight
P = A(1,:)'*Sigma(1,1)*A(1,:);
if rcond(P) == 0
   P = 1.e10*eye(size(A,2));
else 
   P = inv(P);
end 
P
% Initial solution
x = pinv(A(1,:)'*Sigma(1,1)*A(1,:))*A(1,:)'*Sigma(1,1)*b(1)%;

for i = 1:size(b,1)
   K = P*A(i,:)'*inv(A(i,:)*P*A(i,:)'+Sigma(i,i))%;
   P = (eye(size(A,2))-K*A(i,:))*P%;
   x = x+K*(b(i)-A(i,:)*x)%;
%   fprintf('\nSolution:\n');
%   for j = 1:size(A,2)
%      fprintf('  x(%2g) = %6.3f\n',j,x(j));
%   end
end

break

dof = size(b,1)-size(A,2);
if dof ~= 0
   P = (norm(b-A*x))^2*P/dof;
else
   P = (norm(b-A*x))^2*pinv(A'*Sigma*A);  
end      
fprintf('\nFinal Covariance matrix:\n');
for j = 1:size(A,2)
   for k = 1:size(A,2)
      fprintf('%12.7f',P(j,k));
   end
   fprintf('\n');
end
fprintf('\nTrace of Covariance matrix: %12.7f\n',trace(P));
%%%%%%%%%%%%%%%%% end rls.m %%%%%%%%%%%%%%%%%%%




