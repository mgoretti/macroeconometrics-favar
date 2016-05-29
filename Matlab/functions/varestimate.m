function [Y,X,Y_start,Y_hat,PI,e,cons]=varestimate(Data,p,constant)
% The function [Y,X,Y_start,Y_hat,PI,e,cons]=varestimate(Data,p,constant)
% makes use of the SUR representation to estmiate the coefficients,
% the fitted values and the residuals of a VAR(p).
%
% INPUTS: 
% Data - t*n dataset at hand
% p - number of lags 
% constant - ==1 (==0) includes (excludes) the constant in the VAR
%
% OUTPUTS:
% Y - T*n matrix of the SUR representation
% X - T*(n*p) matrix of the SUR representation
% Y_start - p*n first p observations for each variable
% Y_hat - T*n fitted values
% PI - n*(n*p) matrix of coefficients [A1 A2 ... Ap]
% e - T*n estimated residuals
% cons - n*1 constant in the model (set =0 if constant==0) 

% STEP 1: Generate Matrixies for SUR Representation
Y=(Data(p+1:end,:));          % lose first p lags 
X = lagmatrix(Data,1:p);      % generate lagged values 
X(1:p,:)=[];                  % Lose one osb for each lag
 
if constant==1;
  X=[ones(size(X,1),1) X];    % add the constant if needed
end
  
% Isolate the p inital elements of Data
Y_start=(Data(1:p,:)); 

% STEP 2: estimate the model 
MLE=((X'*X)\(X'*Y))'; 
if constant ==1;
cons = MLE(:,1);                % Isolate the constant
PI = MLE(:,2:end);              % Isolate the coefficients 
elseif constant ==0;
cons=zeros(size(MLE,1),1);      % set the constant = 0
PI = MLE;
end

% STEP 3: find fitted values and residuals:
Y_hat=X*MLE';
e = Y-Y_hat;

end