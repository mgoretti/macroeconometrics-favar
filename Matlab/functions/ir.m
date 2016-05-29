function [C,C_plot,Omega_ir,e_chol]=ir(Data,p,constant,hor,varargin)
% The function [C,C_plot,Omega_ir,e_chol]=ir(Data,p,constant,hor,varargin)
% computes structural impulse response functions. If 'cholimpact' is
% specified, it implements short-run restrictions susin a cholesky
% decomposition, otherwise it imposes no restrictions.
%
% INPUT:
% Data     - t*n dataset at hand
% p        - number of lags 
% constant - ==1 (==0) includes (excludes) the constant in the VAR
% hor      - time horizon for the impulse response 
% varagin  - string 'cholimpact' for cholesky on impact (default); 'none' no
%            resrt. In case the string is misspecified or missing it is set
%            to default.
%
% OUTPUT:
% C        - n*n*hor array of impulse responses
% C_plot   - hor*(n*n) reshaped version of C for easier plotting
% Omega_ir - var cov matrix
% e_irChol - structural errors computed as e_chol = S^-1 e_ir

% Step 1: Estimate the model.
[~,~,~,~,PI_ir,e_ir,~]=varestimate(Data,p,constant);

[t,n]=size(Data);        % Number of observations 
T=t-p;                   % Data after losing lags

% STEP 2: Compute the var cov matrix Omega
Omega_ir=(e_ir'*e_ir)./T;

% STEP 3: Insert appropriate restrictions:
if isempty(varargin)==0 & strcmp(varargin{1},'cholimpact') 
    % Cholesky on impact restriction
    H=chol(Omega_ir,'lower'); 
elseif isempty(varargin)==0 & strcmp(varargin{1},'none')
    % No restriction
    H=eye(n);
else
    % Cholesky on impact (default)
    H=chol(Omega_ir,'lower'); 
    warning('Misspecified or missing restriction. Set by default to cholimpact.')
end

% STEP 4:
% Find the companion form matrix
A=companion(PI_ir);

% STEP 5: Create a 3d array to store the powers of A in each layer
wald=zeros(size(A,1),size(A,2),hor);
for i=1:hor 
wald(:,:,i)=A^(i-1);
end

% STEP 6: Select the n*n upper left elements to get the Cj of the wald representation
C=wald(1:n,1:n,:);

% STEP 7: Multiply times S 
for i=1:hor
    C(:,:,i)=C(:,:,i)*H;   
end

% STEP 8: Rearrange the ir in a bidimenstional matrix easier to plot
C_plot=reshape(permute(C, [3 2 1]),[],n*n);

% STEP 9: Compute structural errors:
e_chol=(H\e_ir')';
