% L_DETAIL Step by step explainations of a Matlab implementation of the Lambda
%                   method.  We follow the notation introduced in Strang and Borre, 
%                   pages 495--499

% Written by Kai Borre
% March 23, 2002

fprintf('\n')
echo on
% Test example
I_hat = [5.45;3.1;2.97];
%I_hat =[0; 1];
n = size(I_hat,1);
Sigma = [6.29 5.978 .544; 5.978 6.292 2.34; .544 2.34 6.288];
%Sigma = [53.40 38.40;38.40 28.00];
echo off

% Given float ambiguities
fprintf('\nFloat ambiguities I_hat')
for i = 1:n
    fprintf('\n%12.3f', I_hat(i,1))
end

% Original covariance matrix for I_hat
fprintf('\n\nCovariance matrix for I_hat Sigma') 
for i = 1:n
    fprintf('\n')
    for j = 1:n
        fprintf('%12.3f', Sigma(i,j))
    end
end

fprintf('\n\n')
echo on
% Shift of I_hat in order to secure that -1 < I_hat <= +1
echo off

shifts = I_hat-rem(I_hat,1);
fprintf('\n\nInteger shifts of I_hat')
for i = 1:n
    fprintf('\n%12.0f', shifts(i,1))
end

% Remainders of I_hat
I_hat = rem(I_hat,1);
fprintf('\n\nRemainders of I_hat')
for i = 1:n
    fprintf('\n%12.3f', I_hat(i,1))
end

% L^T*D*L factorization of Sigma
[L,D] = ldldecom(Sigma);
fprintf('\n\nSigma is factorized into L^T*D*L')
fprintf('\n\nThe lower triangular L')
for i = 1:n
    fprintf('\n')
    for j = 1:n
        fprintf('%12.3f',L(i,j))
    end
end
fprintf('\n\nThe diagonal matrix D\n')
for i = 1:n
    fprintf('%12.3f',D(i))
end

% Computing the size of the search volume
chi2 = chistart(D,L,I_hat,1);
fprintf('\n\nchi^2 = %5.3f\n',chi2)

% Doing the decorrelation
[Sigma_t,Z,L_t,D_t] = decorrel(Sigma,I_hat);
fprintf('\nInteger transformation matrix Z')
for i = 1:n
    fprintf('\n')
    for j = 1:n
        fprintf('%12.0f',Z(i,j))
    end
end

% I_hat transformed: z = Z'*I_hat
fprintf('\n\nTransformed, shifted ambiguities z = Z^T*I_hat')
z = Z'*I_hat;
for i = 1:n
    fprintf('\n%12.3f',z(i))
end

% The transformed, decorrelated covariance matrix : Sigma_t = Z^T*Sigma*Z
fprintf('\n\nThe transformed, decorrelated covariance matrix Sigma_t = Z^T*Sigma*Z.')
fprintf('\n\nSigma_t = Z^T*Sigma*Z')
for i = 1:n
    fprintf('\n')
    for j = 1:n
        fprintf('%12.3f',Sigma_t(i,j))
    end
end
fprintf('\nNote the deminished off-diagonal terms!')
fprintf('\n\nThe lower triangular L_t used in the search')
for i = 1:n
    fprintf('\n')
    for j = 1:n
        fprintf('%12.3f',L_t(i,j))
    end
end
fprintf('\n\nThe diagonal matrix D_t used in the seach\n')
for i = 1:n
    fprintf('%12.3f',D_t(i))
end

fprintf('\n\n')
echo on
% Determining the size of the search volume for the transformed L_t and D_t
echo off

chi2 = chistart(D_t,L_t,z,1); 
fprintf('\n\nchi^2 for the transformed problem  =  %5.3f\n',chi2)

fprintf('\nThe search domain is defined as')
fprintf('\n          (I-I_hat)^T*(Z^T*Sigma*Z)*(I-I_hat) < chi^2,  for I integer')

[I_bar,sqnorm,ierr] = lsearch(z,L_t,D_t,chi2,1);
I_bar = (I_bar' * inv(Z))'+shifts;
fprintf('\n\nFixed ambiguities I_bar')
for i = 1:n
    fprintf('\n%12.0f', I_bar(i,1))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% end l_detail.m  %%%%%%%%%%%%
