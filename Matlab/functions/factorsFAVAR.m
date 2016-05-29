function [F] = factorsFAVAR(Data,nf)
% The function F = factorsFAVAR(Data,nf) selects a number of principal
% components of a data set. The function selects the nf first
% eigenvectors, associated to the eigenvlues in descending order, and
% summarizes the information contained in Data in its nf principal
% components.
%
% INPUTS:
% Data - T*N data sample
% nf   - the number of factors to be selected
%
% OUTPUTS:
% F    - T*nf factors 

Gamma = cov(Data);                  % Estimate the varcov matrix of the data
[An,~] = eigs(Gamma,nf,'LM');       % Select the first nf biggest eigenvectors
F = (An'*Data')';                   % Compute the factors

end