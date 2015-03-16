function [H, Qaj] = unidecor(example)
%UNIDECOR Computes a decorrelating transformation as described in
%           Liu et al. (1999) A new approach to GPS ambiguity 
%           decorrelation, Journal of Geodesy, vol 73: 478-490
% Qa......: ambiguity covarinace matrix
% dim.....: dimension of ambiguity covariance matrix
% H.......:	decorrelating transformation
% Qaj.....:	ambiguity correlation matrix in iteration j
% Hj......: decorrelating transformation in iteration j
% Lj......: permutation matrix
% Kj......: lower diagonal matrix, with one column
%				of coefficients of the L-matrix of the
%				LDL'-decomposition
% Final transformation and decorrelated covariance matrix are
% stored in H and Qaj

% Three examples of ambiguity covariance
% matrices, mentioned in the paper
% specify covariance matrix by using the "example"-variable

%Written by Niels Jonkman
%November 1999

if example == 1
   Qa = [80.0 9.75 27.8;
      9.75 1.25 3.45;
      27.8 3.45 10.0];
elseif example == 2
   Qa = [81.708 34.09 75.934 78.342 94.720 146.448;
      34.090 19.611 31.422 34.344 49.839 64.905;
      75.934 31.422 70.641 72.615 87.351 135.801;
      78.342 34.344 72.615 76.029 94.675 141.869;
      94.720 49.839 87.351 94.675 131.472 177.552;
      146.448 64.905 135.801 141.869 177.552 265.416];
else
   Qa = [35.366 40.078 32.142 43.475 42.154 31.859 36.369 39.035 42.985 41.515 45.805 40.521;
      40.078 45.830 36.807 49.742 47.882 36.161 40.784 44.697 48.435 46.933 52.121 45.521; 
      32.142 36.807 30.163 40.284 39.284 29.866 32.807 36.167 38.885 37.503 41.834 36.410;
      43.475 49.742 40.284 56.096 53.512 39.783 44.399 48.693 52.019 51.001 56.903 49.695;	
      42.154 47.882 39.284 53.512 53.400 39.476 44.710 46.956 51.463 49.903 55.093 48.131;
      31.859 36.161 29.866 39.783 39.476 30.295 33.573 35.344 38.779 37.203 42.310 36.344;	
      36.369 40.784 32.807 44.399 44.710 33.573 41.148 38.815 44.450 42.904 46.782 41.726;
      39.035 44.697 36.167 48.693 46.956 35.344 38.815 44.307 47.326 45.624 50.415 44.099; 		
      42.985 48.435 38.885 52.019 51.463 38.779 44.450 47.326 53.096 50.840 55.776 49.417;
      41.515 46.933 37.503 51.001 49.903 37.203 42.904 45.624 50.840 49.216 53.950 47.831;
      45.805 52.121 41.834 56.903 55.093 42.310 46.782 50.415 55.776 53.950 64.832 52.595;
      40.521 45.521 36.410 49.695 48.131 36.344 41.726 44.099 49.417 47.831 52.595 47.218];
end
t = cputime;

%%for q = 1:1000
   % initialization
   % Hj is intialized as a matrix of ones rather than as
   % a unit matrix, in order to get the iteration going
   dim = length(Qa);
   H = eye(dim);
   Hj = ones(dim);
   Qaj = Qa;
   % iterate the HL-processes until the
   % Hj-matrix equals a unit matrix
   while ~isempty(find(Hj-eye(dim)))
      % initialization
      Hj = eye(dim);
      % HL-process
      for i = 1:1:dim-1
         % initialization
         Lj = eye(dim);
         Kj = eye(dim);
         % determine the index of the largest
         % diagonal element of the vc-matrix
         swapindex = find(diag(Qaj) == min(diag(Qaj(i:dim,i:dim))));
         % build permutation matrix
         Lj(i,i) = 0;
         Lj(swapindex,swapindex) = 0;
         Lj(i,swapindex) = 1;
         Lj(swapindex,i) = 1;
         % permutate covariance matrix
         Qaj = Lj*Qaj*Lj';
         % determine column i of L-matrix in
         % LDL'-decomposition and store in matrix Kj
         Kj(i+1:dim,i) = -Qaj(i+1:dim,i)/Qaj(i,i);
         % apply transformation
         Qaj = Kj*Qaj*Kj';
         % determine "float" transformation matrix
         Hj = Kj*Lj*Hj;
      end;
      % round the "float" transformation matrix
      Hj = round(Hj);
      % determine transformation matrix after iteration j
      H = Hj*H;
      % apply transformation
      Qaj = H*Qa*H';
   end;
%% end
cputime-t
%%%%%%%%%%%%%%%%%%%% end unidecor.m  %%%%%%%%%%%%%