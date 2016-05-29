function [comp]=companion(PI)
% The function [comp]=companion(PI) computes the companion form of a matrix
% of coefficients. 
%
% INPUTS:
% PI   - matrix of coefficents
%
% OUTPUTS :
% comp - n*p square matrix in the companion form

s=size(PI,2);
comp = (eye(s));
comp=[PI; comp];
comp=comp(1:s,:);

end