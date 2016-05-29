function [ULb,FEVULB,IR,FEVm,ULb2,F4]=bootbands(Data,p,constant,iter,hor,conf,varargin)
% The function [ULb,FEVULB,IR,FEV,ULB2]=bootbands(Data,p,constant,iter,hor,conf,varargin)
% computes structural impulse response functions and variance decomposition
% with bootstrap bands.
%
% INPUTS:
% Data     - t*n dataset at hand
% p        - number of lags
% constant - ==1 (==0) includes (excludes) the constant in the VAR 
% iter     - number of bootstrap iteration
% hor      - time horizon for the impulse response
% conf     - confidence level for the bands / eiter scalar or 2-vector
% varagin  - string 'cholimpact' for cholesky on impact (default); 'none' no
%            resrt. In case the string is misspecified or missing it is set
%            to default.
%
% OUTPUTS:
% ULb      - n*n*hor*2 upper and lower bands in tensor format
% FEVULB   - n*n*hor*2 upper and lower 1std bands in tensor format
% IR       - n*n*hor impulse responses (see function below)
% FEVm     - n*n*hor mean Variance decompostion (see function below)
% ULb2     - n*n*hor*2 upper and lower bands if there are 2 conf specified
% F4       - n*n*hor*iter all IR of the bootstrap


% Step 1: Estimate the model.
[~,~,Y_start,~,PI,e,cons]=varestimate(Data,p,constant);

[t,n]=size(Data);        % Number of observations 
T=t-p;                   % Data after losing lags

% Compute IR and Var decomposition:
IR=ir(Data,p,constant,hor,'cholimpact');

% Prepare matrixes for looping
F4=zeros(n,n,hor,iter);
FEV4=zeros(n,n,hor,iter);

for j=1:iter
% STEP 2: generate a new sample with bootstrap
Y_new=zeros(t-p,n);
Y_in= reshape(flipud(Y_start)',1,[]);                  % 1*(n*p) row vector of [y_0 ... y_-3]
  for i=1:T;
    Y_new(i,:)= cons' + Y_in*PI' + e(randi(T),:);  
    Y_in=[Y_new(i,:) Y_in(1:n*(p-1))];
  end
    Data_new=[Y_start ; Y_new];                        % Add the p initial lags to recreate the sample
  
% STEP 3: Do impulse response
if isempty(varargin)==0 && strcmp(varargin{1},'cholimpact') 
    % Cholesky on impact restriction
    IR_new=ir(Data_new,p,constant,hor,'cholimpact');
elseif isempty(varargin)==0 && strcmp(varargin{1},'none')
    % No restriction
    IR_new=ir(Data_new,p,constant,hor,'none');
else
    % Cholesky on impact (default)
    IR_new=ir(Data_new,p,constant,hor,'cholimpact');
    %warning('Misspecified or missing restriction. Set by default to cholimpact.')
end

% STEP 4: Store in a 4d array:
F4(:,:,:,j)=IR_new;
FEV4(:,:,:,j)=vardec(IR_new); % Compute variance decompostion

end
% STEP 5: Find the percentile for the IR bands:
ULb=prctile(F4,[(100-conf(1))/2 (100+conf(1))/2],4);
ULb2=NaN(size(ULb));
if length(conf)>1
    ULb2=prctile(F4,[(100-conf(2))/2 (100+conf(2))/2],4);
end

% STEP 6: Find mean and 1 std for the variance dec.
FEVm=mean(FEV4,4);
%FEVm_plot=reshape(permute(FEVm,[3 2 1]),[],n^2);
FEVstd=std(FEV4,0,4);
FEVULB(:,:,:,1)=FEVm+FEVstd;
FEVULB(:,:,:,2)=FEVm-FEVstd;
%FEVulb=reshape(permute(FEVULB,[3 2 1 4]),[],n^2,2);
end

function [FEV, FEV_plot]=vardec(IR)
% The function FEV=vardec(IR) takes a set of impulse responses and computes
% the variance decomposition along all its horizons:
% 
% INPUTS;
% IR       - n*n*hor tensor of impulse responses whose shocks have UNIT VARIANCE
%
% OUTPUT: 
% FEV      - n*n*hor tensor of explained variance at each horizon
% FEV_plot - hor*(n*n) reshaped version of FEV for easy plotting

% Read matrix size;
n=size(IR,1);

% Cumulate the squared sum of coeffs:
F2=cumsum(IR.^2,3);
% Total repeated for each entry of the tensor
totVAR=repmat(sum(F2,2),1,n); 
% Variance decomposition
FEV=F2./totVAR;
FEV_plot=reshape(permute(FEV,[3 2 1]),[],n*n);
end

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

end
