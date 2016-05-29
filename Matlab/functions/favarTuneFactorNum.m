function [ facCVrel, CV, CVabs, CVrel ] = favarTuneFactorNum(Data, svarData, lags, maxFactors, constant )
% finds the optimal number of factors for the FAVAR by finding the
% lowest test error of the cross-validation
% 
% INPUTS:
% Data: dataset used for PCA
% svarData: dataset of the variables that will be also used as the explained variables
% lags: number of lags in the VAR model
% maxFactors: maximal number of factors to test
% constant: boolean deciding whether to have a constant term in the VAR model
%
% OUTPUTS:
% facCVrel: optimal number of factors
% CV: raw CV RMSE
% CVabs: mean of the RMSE of each variable
% CVrel: ponderated mean of the RMSE of each variable
 
setSeed(1)

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

for f = 1:maxFactors   
    F = factorsFAVAR(Data, f);
    favarData = [svarData F];
    CV(f, :) = CVFavar(favarData, lags, f, constant, idxCV);
end

CVabs = sum(CV, 2);

for i=1:size(CV,1)
    CVrel(i, :) = CV(i, :) ./ mean(CV,1);
end
CVrel = sum(CVrel, 2);
[~,facCVrel] = min(CVrel); % # of factors suggested by CVrel

end

