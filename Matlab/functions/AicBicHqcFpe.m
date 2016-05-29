function [InformationCriterion, aicL, bicL, hqcL, fpeL, lagsCVabs, lagsCVrel, CVrelTr] = AicBicHqcFpe(Data, pmax)
% The function [InformationCriterion, aicL, bicL, hqcL, fpeL] = AicBicHqcScFpe(Data, pmax) 
% gives the AIC, BIC, HQIC and FPE information criteria for lag selection.
%
% Akaike information criterion :
%     AIC = logV + 2d/T 
% Bayesian information criterion
%     BIC = logV + logT*d/T
% Hannan-Quinn information criterion
%     HQC = logV + 2*d/T*log(log(T))
% Final Predicton Error criterion
%     FPE = V * ((T + Np + 1)/(T - Np - 1))^N
%
% where:
% V = det(1/T Sum_{t=1}^T e_t*e_t')
% d : # of parameters 
% T : # of observations
% p : # of lags
%
% INPUTS:
% Data : T*N data sample
% pmax : scalar maximum number of lags to be checked
%
% OUTPUTS:
% InformationCriterion: raw scores of each test for each lag
% aicL : number of lags suggested by AIC
% bicL : number of lags suggested by BIC
% hqcL : number of lags suggested by HQIC
% fpeL : number of lags suggested by FPE
% lagsCVabs: number of lags suggested by CV (non-scaled by mean error for the variable)
% lagsCVrel: number of lags suggested by CV (scaled by mean error for the variable)
% CVrelTr: raw training error of the CV (testing error included in InformationCriterion)


% generate M distributions of the splits for the CVs
setSeed(1);
M = 100;
K = 5;
N = size(Data, 1);
Nk = floor(N/K);

idxCV = NaN (K, Nk, M);
for m = 1:M;
    idx = randperm(N);
    for k = 1:K
        idxCV(k,:, m) = idx(1+(k-1)*Nk:k*Nk);
    end
end

% Read Data:
[T,N] = size(Data);
constant=1;
% Preallocate for looping:
aic = NaN(pmax,1);
bic = NaN(pmax,1);
hqc = NaN(pmax,1);
fpe = NaN(pmax,1);
CVTe = NaN(pmax,size(Data, 2));
CVTr = NaN(pmax,size(Data, 2));
for p = 1:pmax
    [~,~,~,~,~,e]=varestimate(Data,p,constant); % Estimate the model with p lags
    V = det((e'*e)/T);                          % Estimate the varcov matrix 
    d = p*N^2;                                  % Number of parameters
    aic(p) = log(V) + 2*d/T;                    % Estimate AIC 
    bic(p) = log(V) + log(T)*d/T;               % Estimate BIC 
    hqc(p) = log(V) + 2*d/T*log(log(T));        % Estimate HQ 
    fpe(p) = V*((T + N*p + 1)/(T - N*p - 1))^N; % Estimate FPE 
    [CVTe(p, :), CVTr(p, :)] = CVvar(Data, p, constant, idxCV);

    % Check that the varcov is positive definite
    [~,test]=chol(e'*e);
    if test~=0 % if not discard the lag option
        aic(p)=Inf; bic(p)=Inf; hqc(p)=Inf; fpe(p)=Inf;
    end
end
% Matrix with estimated criterions, rows are number of lags
[~,aicL] = min(aic) % # of lags suggested by AIC
[~,bicL] = min(bic) % # of lags suggested by BIC
[~,hqcL] = min(hqc) % # of lags suggested by HQIC
[~,fpeL] = min(fpe) % # of lags suggested by FPE

% CV
CVTe

% absolute error (sum of the RMSE of each col)
CVabs = sum(CVTe, 2);
[~,lagsCVabs] = min(CVabs) % # of lags suggested by CVabs


% relative error (sum of the errors scaled by the mean error of the col)
for i=1:size(CVTr,1)
    CVrelTr(i, :) = CVTr(i, :) ./ mean(CVTr,1);
end
CVrelTr = sum(CVrelTr, 2);

for i=1:size(CVTe,1)
    CVrelTe(i, :) = CVTe(i, :) ./ mean(CVTe,1);
end
CVrelTe = sum(CVrelTe, 2);
[~,lagsCVrel] = min(CVrelTe) % # of lags suggested by CVrel

InformationCriterion = [aic bic hqc fpe CVabs CVrelTe]

end
