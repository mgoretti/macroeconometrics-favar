function  screeplot(Data, nf)
% The function screeplot(Data, nf) draws a plot of the percentage of the
% total amount of variability of the original data set explained by a 
% certain number of factors. Two lines are shown, one is the cumulative 
% percent explained and the other is the marginal gain of including the ith
% factor.
%
% INPUTS:
% Data - T*N data sample
% nf   - maximum number of factors on the scree plot

[~, ~, C] = pca(Data);

varTot = sum(C);
p = 100 * C/varTot;
p = p(1:nf,:);
pe = cumsum(p);

line(1:1:nf, pe(1:1:nf),'marker', '.', 'color', 'b', 'markerfacecolor', 'g','LineWidth',1.5);
line(1:1:nf, p(1:1:nf),'marker', '.', 'color', 'k', 'markerfacecolor', 'g','LineWidth',1.5);
xlabel('Number of factors');
ylabel('Percent explained' );
legend( {'Cumulative', 'Individual'}, 'location', 'northwest' );
title( 'Scree plot');

end

