%------------------------------------------------------------------------%
%------------------------------------------------------------------------%
%               Macroeconometrics - Spring 2016
%
%       Macroeconomic effects of oil price shocks : a FAVAR approach.
%
%                       May 28, 2016
%
%                       Marco Goretti
%                     Frederic Martenet
%                       Simon Tieche
%------------------------------------------------------------------------%
%------------------------------------------------------------------------%

% This file is structured as follows :
%    A) Data Preparation
%    B) SVAR
%    C) FAVAR
%    D) Robustness Checks

% Necessary functions :
%    - AicBibHqcFpe.m
%    - bootbands.m
%    - bw_trend.m
%    - calendar_make.m
%    - companion.m
%    - CVFavar.m
%    - CVvar.m
%    - factorsFAVAR.m
%    - favarTuneFActorsNum.m
%    - formatGraph.m
%    - ir.m
%    - onatski.m
%    - saveGraph.m
%    - screeplot.m
%    - setSeed.m
%    - simontable.m
%    - varestimate.m 

%------------------------------------------------------------------------%

clc; clear all; close all;

% not include in the final code
path = '/Users/Marco/Dropbox/Macroeconometrics/Projet/';
% path = '/Users/FredericMartenet/Desktop/Dropbox/FRED/Documents_dropbox/HEC/_MScE/_4_2/Macroeconometrics/Macroeconometrics/Projet/' 
%path = 'C:\Users\Simon\Dropbox\Macroeconometrics\Projet\';



%path = '... path of your directroy here ...'

cd (path)                   % Set working directory
addpath(path,'Functions')   % Set folder with functions

%%
%------------------------------------------------------------------------%
%   A) Data Preparation
%------------------------------------------------------------------------%

%------ Import Data ------%

load('Data/DATA');      
% data_factors : Data set containing 106 variables describing the US
%                economy, from Stock & Watson (2915). Used to estimate 
%                the factors for the FAVAR.
% FF           : Effective Fed Funds Rate, percent.
% GDP          : Real GDP, 3 decimals, billions of Chained 2000 Dollars.
% oil          : Oil price index, index 1982=100.
% PCEPI        : Consumer price index, index 2009=100.
% TOTEMP       : Total nonfarm employment, thousands of Persons.
% All variables are quarterly.

%------ Descriptive Evidence ------%

[dnobs,time,~] = calendar_make([1959 1],[2014 3],4);   % time vector for plots

% Producer Price Index by Commodity for Fuels and Related Products and Power: Crude Petroleum 
% Source    :   US. Bureau of Labor Statistics
% Units     :   Index 1982=100
figure(1)
plot(time,oil,'Color','black','LineWidth',2);
xlim([1959 time(dnobs)]);grid on;
hxlabel = 'Time';
hylabel = 'Oil price,Index 1982=100';
saveGraph(strcat(path,'Paper/Figures/Descriptive_oil'), hxlabel,  hylabel, 1.4);

% Real GDP, 3 decimal
% Source    :   US. Bureau of Economic Analysis
% Units     :   Billions of Chained 2000 Dollars
figure(2)
plot(time,GDP,'Color','black','LineWidth',2);
xlim([1959 time(dnobs)]);grid on;
hxlabel = 'Time';
hylabel = 'Real GDP, 3 decimals';
saveGraph(strcat(path,'Paper/Figures/Descriptive_GDP'), hxlabel,  hylabel, 1.4);

% All Employees: Total Nonfarm Payrolls
% Source    :   US. Bureau of Labor Statistics
% Units     :   Thousands of Persons
figure(3)
plot(time,TOTEMP,'Color','black','LineWidth',2);
xlim([1959 time(dnobs)]);grid on;
hxlabel = 'Time';
hylabel = 'Total Employment, Thousands of Persons';
saveGraph(strcat(path,'Paper/Figures/Descriptive_TOTEMP'), hxlabel,  hylabel, 1.4);

% Personal Consumption Expenditures: Chain-type Price Index
% Source    :   US. Bureau of Economic Analysis
% Units     :   Index 2009=100
figure(4)
plot(time,PCEPI,'Color','black','LineWidth',2);
xlim([1959 time(dnobs)]);grid on;
hxlabel = 'Time';
hylabel = 'Consumer Price Index, Index 2009=100';
saveGraph(strcat(path,'Paper/Figures/Descriptive_PCEPI'), hxlabel,  hylabel, 1.4);

% Effective Federal Funds Rate
% Source    :   Board of Governors of the Federal Reserve System (US)
% Units     :   Percent 
figure(5)
plot(time,FF,'Color','black','LineWidth',2);
xlim([1959 time(dnobs)]);grid on;
hxlabel = 'Time';
hylabel = 'Fed Funds Rate, percent';
saveGraph(strcat(path,'Paper/Figures/Descriptive_FF'), hxlabel,  hylabel, 1.4);

close all;
clearvars dnobs time hxlabel hylabel

%------ Data Transformation ------%

% Oil price : first differences of log
oil = diff(log(oil));
oil(1,:) = [];

% GDP : first differences of log
gdp = diff(log(GDP));
gdp(1,:)=[];

% Total employment : differences of log
emp = diff(log(TOTEMP));
emp(1,:) = [];

% PCEPI : second differences of log
pcepi = diff(diff(log(PCEPI))); 

% Fed Funds Rate : first differences
ff = diff(FF);
ff(1,:) = [];

% Augmented DF tests for unit root
[oil_yn ,oil_p, oil_tstat]       = adftest(oil);
[gdp_yn ,gdp_p, gdp_tstat]       = adftest(gdp);
[emp_yn ,emp_p, emp_tstat]       = adftest(emp);
[ff_yn ,ff_p, ff_tstat]          = adftest(ff);
[pcepi_yn ,pcepi_p, pcepi_tstat] = adftest(pcepi);
% Table with results
DF = [oil_tstat gdp_tstat emp_tstat pcepi_tstat ff_tstat;... % t statistics
      oil_p     gdp_p     emp_p     pcepi_p     ff_p;...     % p-value
      oil_yn    gdp_yn    emp_yn    pcepi_yn    ff_yn];      % rejected Y or N
name = 'Paper/Tables/DF.tex'; % for Mac
%name = 'Paper\Tables\DF.tex'; % for PC
caption = 'Augmented Dickey-Fuller tests for unit root';
numbercolumn = 'lccccc';
headers = ' & Oil price & Real GDP & Employment & Inflation & Fed Funds rate'; % Warning ! Must have the same length
namecolumn = ['Test statistic                 '; 'P-value                        '; 'Rejection of the unit-root null'];
note = 'This table reports the final prediction error (FPE), Akaike''s information criterion (AIC), Swarz''s Bayesian information criterion (BIC), and the Hannan and Quinn information criterion (HQIC) lag-order selection statistics for a series of vector autoregressions of order 1 through a maximum of 12 lags. ';
tablewidth = '1';
label = 'tab:DF';
simontable(DF, name, label, caption, tablewidth, numbercolumn, namecolumn, headers, note)

clearvars DF oil_tstat gdp_tstat emp_tstat pcepi_tstat ff_tstat... 
      oil_p     gdp_p     emp_p     pcepi_p     ff_p...   
      oil_yn    gdp_yn    emp_yn    pcepi_yn    ff_yn...
      caption headers numbercolumn note name tablewidth label namecolumn

%------ Data set containing all variables ------%

DataSVAR=[oil gdp emp pcepi ff];
nt = size(DataSVAR,1);           % number of time periods
ns = size(DataSVAR,2);           % number of series

%------ Low frequency trends ------%

% Removing low frequency trends using a Bi-Weight trend
bw_bw = 100; % Bi-Weight Parameter for local demeaning
for is = 1:ns; 
    tmp = bw_trend(DataSVAR(:,is),bw_bw);
   	DataSVAR_trend(:,is)= tmp;
   	DataSVAR(:,is) = DataSVAR(:,is) - DataSVAR_trend(:,is); 
end;

%------ Standardization ------%
xmean = nanmean(DataSVAR)';                                        % mean (ignoring NaN)
mult = sqrt((sum(~isnan(DataSVAR))-1)./sum(~isnan(DataSVAR)));     % num of non-NaN entries for each series
xstd = (nanstd(DataSVAR).*mult)';                                  % std (ignoring NaN)
DATA_SVAR = (DataSVAR - repmat(xmean',nt,1))./repmat(xstd',nt,1);  % standardized data

clearvars tmp DataSVAR DataSVAR_trend xmean xstd mult is bw_bw

% Final Data Set : DATA_SVAR
% It contains the 5 variables of interest ready for the SVAR estimation.


%%
%------------------------------------------------------------------------%
%   B) SVAR
%------------------------------------------------------------------------%

%------ Lag Selection ------%

[InformationCriterion, aicL, bicL, hqcL, fpeL, CVabs, CVrel, CVrelTr] = AicBicHqcFpe(DATA_SVAR, 12);

%------ Cross Validation results ---------%
CVrelTe = InformationCriterion(:, end);
plot(1:length(CVrelTr), CVrelTr, 'color', [77 45 115]/256,'LineWidth',2);
hold on; 
plot(1:length(CVrelTr), CVrelTe, 'color', [17 152 152]/256,'LineWidth',2);
grid on;
legend('Train','Test', 'Location','northeast')
xlim([1 12]);

hx = xlabel('lags');
hy = ylabel('mean RMSE');

set(gca,'fontsize',14,'fontname','Helvetica','box','off','tickdir','out','ticklength',[.02 .02],'xcolor',0.5*[1 1 1],'ycolor',0.5*[1 1 1]);
set([hx; hy],'fontsize',12,'fontname','avantgarde','color',[.3 .3 .3]);
%grid on;

hold off;
w = 7; h = 5;
set(gcf, 'PaperPosition', [0 0 w h]); %Position plot at left hand corner with width w and height h.
set(gcf, 'PaperSize', [w h]); %Set the paper to have width w and height h.
saveas(gcf, strcat(path,'Paper/Figures/CVlags'), 'pdf') %Save figure

close all;

% Table with results
name = 'Paper/Tables/IC.tex'; % for Mac
%name = 'Paper\Tables\IC.tex'; % for PC
caption = 'Information Criteria';
numbercolumn = 'lcccccc';
namecolumn = ['1 '; '2 '; '3 '; '4 '; '5 '; '6 '; '7 ';'8 '; '9 ';'10';'11';'12']; % Warning ! Must have the same length
headers = 'lag & AIC & BIC & HQIC & FPE & CVabs & Cvrel';
note = 'This table reports the final prediction error (FPE), Akaike''s information criterion (AIC), Swarz''s Bayesian information criterion (BIC), the Hannan and Quinn information criterion (HQIC), Cross validation criterion (absolute CVabs and relative CVrel) lag-order selection statistics for a series of vector autoregressions of order 1 through a maximum of 12 lags. ';
tablewidth = '0.55';
label = 'tab:IC';
simontable(InformationCriterion, name, label, caption, tablewidth, numbercolumn, namecolumn, headers, note)

clearvars InformationCriterion aicL bicL hqcL fpeL CVabs CVrelTr CVrelTe...
    caption headers numbercolumn note tablewidth label namecolumn name

%------ SVAR estimation ------%

% Parameters        
p = CVrel;                 % number of lags
constant = 1;              % include the constant
hor = 20;                  % horizon
iter = 500;                % number of iterations
conf = [90 60];            % level for confidence bands

% Estimation of the SVAR
Cf_SVAR = ir(DATA_SVAR,p,constant,hor,'cholimpact');
[ULbf_SVAR,~,~,~,ULb2f_SVAR] = bootbands(DATA_SVAR,p,constant,iter,hor,conf,'cholimpact');

% Plot resutls
irfs_SVAR = permute(Cf_SVAR,[3 1 2]);
irfs_SVAR = irfs_SVAR(:,:,1);
conf_SVAR1 = permute(ULbf_SVAR,[3 1 2 4]);
conf_SVAR_down1 = conf_SVAR1(:,:,1,1);
conf_SVAR_up1 = conf_SVAR1(:,:,1,2);
conf_SVAR2 = permute(ULb2f_SVAR,[3 1 2 4]);
conf_SVAR_down2 = conf_SVAR2(:,:,1,1);
conf_SVAR_up2 = conf_SVAR2(:,:,1,2);
plot_titles = { 'Oil Price', 'Real GDP' , 'Total Employment', 'Inflation'};

figure(6);
for ii = 1:4
subplot(2,2,ii); 
% IRFs
plot(irfs_SVAR(:,ii),'black','LineWidth',2);title(plot_titles(:,ii));
formatGraph('Quarter', '');
% Horizontal line at 0
hold on;plot(xlim, [0 0], 'color', 0.5*[1 1 1],'LineWidth',1.1);                     
% 90% confidence interval
hold on;plot(conf_SVAR_down1(:,ii),'black','LineStyle',':');hold on;plot(conf_SVAR_up1(:,ii),'black','LineStyle',':');
% 60% confidence interval
hold on;plot(conf_SVAR_down2(:,ii),'black','LineStyle','-');hold on;plot(conf_SVAR_up2(:,ii),'black','LineStyle','-');
end
saveGraph(strcat(path,'Paper/Figures/SVAR_irf'), 'Quarter',  '', 1.4);
close all;

clearvars plot_titles constant hor iter conf...
    conf_SVAR1 conf_SVAR2 conf_SVAR_down1 conf_SVAR_down2...
    conf_SVAR_up1 conf_SVAR_up2 ii ULbf_SVAR ULb2f_SVAR...
    constant hor iter conf

%%
%------------------------------------------------------------------------%
%   C) FAVAR
%------------------------------------------------------------------------%

%------ Number of factors ------%

% Screeplot 
screeplot(data_factors, 60);
hxlabel = 'Number of factors';
hylabel = 'Percent explained';

hx = xlabel(hxlabel);
hy = ylabel(hylabel);
grid on;
set(gca,'fontsize',14,'fontname','Helvetica','box','off','tickdir','out','ticklength',[.02 .02],'xcolor',0.5*[1 1 1],'ycolor',0.5*[1 1 1]);
set([hx; hy],'fontsize',12,'fontname','avantgarde','color',[.3 .3 .3]);
%grid on;

hold off;
w = 5*1.4; h = 5;
set(gcf, 'PaperPosition', [0 0 w h]); %Position plot at left hand corner with width w and height h.
set(gcf, 'PaperSize', [w h]); %Set the paper to have width w and height h.
saveas(gcf, strcat(path,'Paper/Figures/screeplot'), 'pdf') %Save figure


%saveGraph(strcat(path,'Paper/Figures/screeplot'), hxlabel, hylabel, 1.4);
close all;
clearvars hxlabel hylabel

% Cross Validation
constant = 1;
[numFt, FAVAR_CV_RMSE, CVabs, CVrel] = favarTuneFactorNum(data_factors, DATA_SVAR(:, :), p, 12, constant)

% Onatski(2010)
numFtOnatski = onatski(data_factors,12);

numFt == numFtOnatski     % Both indicate 2 factors

%------ Factors estimation ------%

F = factorsFAVAR(data_factors,numFt);

%------ FAVAR estimation ------%

% Data
DATA_FAVAR = [DATA_SVAR F];

% Parameters
%p = CVrel;                  % number of lags
constant = 1;              % include the constant
hor = 20;                  % horizon
iter = 500;                % number of iterations
conf = [90 60];            % level for confidence bands

% Estimation of the FAVAR
Cf_FAVAR = ir(DATA_FAVAR,p,constant,hor,'cholimpact');
[ULbf_FAVAR,~,~,~,ULb2f_FAVAR] = bootbands(DATA_FAVAR,p,constant,iter,hor,conf,'cholimpact');

% Plot resutls
irfs_FAVAR = permute(Cf_FAVAR,[3 1 2]);
irfs_FAVAR = irfs_FAVAR(:,:,1);
conf_FAVAR1 = permute(ULbf_FAVAR,[3 1 2 4]);
conf_FAVAR_down1 = conf_FAVAR1(:,:,1,1);
conf_FAVAR_up1 = conf_FAVAR1(:,:,1,2);
conf_FAVAR2 = permute(ULb2f_FAVAR,[3 1 2 4]);
conf_FAVAR_down2 = conf_FAVAR2(:,:,1,1);
conf_FAVAR_up2 = conf_FAVAR2(:,:,1,2);
plot_titles = { 'Oil Price', 'Real GDP' , 'Total Employment', 'Inflation'};

figure(7);
for ii = 1:4
subplot(2,2,ii); 
% IRFs of the SVAR
plot(irfs_SVAR(:,ii),'red','LineWidth',1.5,'LineStyle','-');hold on; 
% IRFs of the FAVAR
plot(irfs_FAVAR(:,ii),'black','LineWidth',2);title(plot_titles(:,ii));  
formatGraph('Quarter', '');
% Horizontal line at 0
hold on;plot(xlim, [0 0], 'color', 0.5*[1 1 1],'LineWidth',1.1);                     
% 90% confidence interval
hold on;plot(conf_FAVAR_down1(:,ii),'black','LineStyle',':');hold on;plot(conf_FAVAR_up1(:,ii),'black','LineStyle',':')
% 60% confidence interval
hold on;plot(conf_FAVAR_down2(:,ii),'black','LineStyle','-');hold on;plot(conf_FAVAR_up2(:,ii),'black','LineStyle','-')
end
saveGraph(strcat(path,'Paper/Figures/FAVAR_irf'), 'Quarter',  '', 1.4);
close all;

clearvars plot_titles constant hor iter conf...
    irfs_FAVAR conf_FAVAR1 conf_FAVAR2 conf_FAVAR_down1 conf_FAVAR_down2...
    conf_FAVAR_up1 conf_FAVAR_up2 ii numFtOnatski  constant hor iter conf

%%
%------------------------------------------------------------------------%
%   D) Robustness Checks
%------------------------------------------------------------------------%

%------------------------------------------------------------------------%
%------ Restict sample: without zero lower bound (until 2008:Q3) ------%
%------------------------------------------------------------------------%
%------ Sample Selection ------%
% SVAR
DataSVAR = [oil gdp emp pcepi ff];
DataSVAR = DataSVAR(1:199,:);
nt = size(DataSVAR,1);         
ns = size(DataSVAR,2);           
bw_bw = 100; 
for is = 1:ns; 
    tmp = bw_trend(DataSVAR(:,is),bw_bw);
   	DataSVAR_trend_robust1(:,is)= tmp;
   	DataSVAR(:,is) = DataSVAR(:,is) - DataSVAR_trend_robust1(:,is); 
end;
xmean = nanmean(DataSVAR)';                                       
mult = sqrt((sum(~isnan(DataSVAR))-1)./sum(~isnan(DataSVAR)));     
xstd = (nanstd(DataSVAR).*mult)';                                  
DATA_SVAR_robust1 = (DataSVAR - repmat(xmean',nt,1))./repmat(xstd',nt,1);  
% FAVAR
DATA_FACTORS = data_factors(1:199,:);
F = factorsFAVAR(DATA_FACTORS,numFt);
DATA_FAVAR_robust1 = [DATA_SVAR_robust1 F];
%------ Estimation ------%
% SVAR
%p = hqcL;                  
constant = 1;             
hor = 20;                
iter = 500;                
conf = [90 60];            
Cf_SVAR_robust1 = ir(DATA_SVAR_robust1,p,constant,hor,'cholimpact');
[ULbf_SVAR_robust1,~,~,~,ULb2f_SVAR_robust1] = bootbands(DATA_SVAR_robust1,p,constant,iter,hor,conf,'cholimpact');
irfs_SVAR = permute(Cf_SVAR_robust1,[3 1 2]);
irfs_SVAR = irfs_SVAR(:,:,1);
% FAVAR
Cf_FAVAR_robust1 = ir(DATA_FAVAR_robust1,p,constant,hor,'cholimpact');
[ULbf_FAVAR_robust1,~,~,~,ULb2f_FAVAR_robust1] = bootbands(DATA_FAVAR_robust1,p,constant,iter,hor,conf,'cholimpact');
irfs_FAVAR = permute(Cf_FAVAR_robust1,[3 1 2]);
irfs_FAVAR = irfs_FAVAR(:,:,1);
conf_FAVAR1 = permute(ULbf_FAVAR_robust1,[3 1 2 4]);
conf_FAVAR_down1 = conf_FAVAR1(:,:,1,1);
conf_FAVAR_up1 = conf_FAVAR1(:,:,1,2);
conf_FAVAR2 = permute(ULb2f_FAVAR_robust1,[3 1 2 4]);
conf_FAVAR_down2 = conf_FAVAR2(:,:,1,1);
conf_FAVAR_up2 = conf_FAVAR2(:,:,1,2);
plot_titles = { 'Oil Price', 'Real GDP' , 'Total Employment', 'Inflation'};
figure(8);
for ii = 1:4
subplot(2,2,ii); 
plot(irfs_SVAR(:,ii),'red','LineWidth',1.5,'LineStyle','-');hold on;      
plot(irfs_FAVAR(:,ii),'black','LineWidth',2);title(plot_titles(:,ii)); 
formatGraph('Quarter', '');
hold on;plot(xlim, [0 0], 'color', 0.5*[1 1 1],'LineWidth',1.1);                     
hold on;plot(conf_FAVAR_down1(:,ii),'black','LineStyle',':');hold on;plot(conf_FAVAR_up1(:,ii),'black','LineStyle',':')
hold on;plot(conf_FAVAR_down2(:,ii),'black','LineStyle','-');hold on;plot(conf_FAVAR_up2(:,ii),'black','LineStyle','-')
end
saveGraph(strcat(path,'Paper/Figures/Robustness_ZLB'), 'Quarter',  '', 1.4);

%%
%------------------------------------------------------------------------%
%------ Restict sample: before/after great moderation ------%
%------------------------------------------------------------------------%
%------ BEFORE the great moderation : 1959-1984 ------%
%------ Sample Selection ------%
% SVAR
DataSVAR = [oil gdp emp pcepi ff];
DataSVAR = DataSVAR(1:104,:);
nt = size(DataSVAR,1);         
ns = size(DataSVAR,2);           
bw_bw = 100; 
for is = 1:ns; 
    tmp = bw_trend(DataSVAR(:,is),bw_bw);
   	DataSVAR_trend_robust2a(:,is)= tmp;
   	DataSVAR(:,is) = DataSVAR(:,is) - DataSVAR_trend_robust2a(:,is); 
end;
xmean = nanmean(DataSVAR)';                                       
mult = sqrt((sum(~isnan(DataSVAR))-1)./sum(~isnan(DataSVAR)));     
xstd = (nanstd(DataSVAR).*mult)';                                  
DATA_SVAR_robust2a = (DataSVAR - repmat(xmean',nt,1))./repmat(xstd',nt,1);  
% FAVAR
DATA_FACTORS = data_factors(1:104,:);
F = factorsFAVAR(DATA_FACTORS,numFt);
DATA_FAVAR_robust2a = [DATA_SVAR_robust2a F];
%------ Estimation ------%
% SVAR
%p = hqcL;                  
constant = 1;             
hor = 20;                
iter = 500;                
conf = [90 60];            
Cf_SVAR_robust2a = ir(DATA_SVAR_robust2a,p,constant,hor,'cholimpact');
irfs_SVAR = permute(Cf_SVAR_robust2a,[3 1 2]);
irfs_SVAR = irfs_SVAR(:,:,1);
% FAVAR
Cf_FAVAR_robust2a = ir(DATA_FAVAR_robust2a,p,constant,hor,'cholimpact');
[ULbf_FAVAR_robust2a,~,~,~,ULb2f_FAVAR_robust2a] = bootbands(DATA_FAVAR_robust2a,p,constant,iter,hor,conf,'cholimpact');
irfs_FAVAR = permute(Cf_FAVAR_robust2a,[3 1 2]);
irfs_FAVAR = irfs_FAVAR(:,:,1);
conf_FAVAR1 = permute(ULbf_FAVAR_robust2a,[3 1 2 4]);
conf_FAVAR_down1 = conf_FAVAR1(:,:,1,1);
conf_FAVAR_up1 = conf_FAVAR1(:,:,1,2);
conf_FAVAR2 = permute(ULb2f_FAVAR_robust2a,[3 1 2 4]);
conf_FAVAR_down2 = conf_FAVAR2(:,:,1,1);
conf_FAVAR_up2 = conf_FAVAR2(:,:,1,2);
plot_titles = { 'Oil Price', 'Real GDP' , 'Total Employment', 'Inflation'};
figure(9);
for ii = 1:4
subplot(2,2,ii); 
plot(irfs_SVAR(:,ii),'red','LineWidth',1.5,'LineStyle','-');hold on;    
plot(irfs_FAVAR(:,ii),'black','LineWidth',2);title(plot_titles(:,ii));  
formatGraph('Quarter', '');
hold on;plot(xlim, [0 0], 'color', 0.5*[1 1 1],'LineWidth',1.1);                     
hold on;plot(conf_FAVAR_down1(:,ii),'black','LineStyle',':');hold on;plot(conf_FAVAR_up1(:,ii),'black','LineStyle',':')
hold on;plot(conf_FAVAR_down2(:,ii),'black','LineStyle','-');hold on;plot(conf_FAVAR_up2(:,ii),'black','LineStyle','-')
end
saveGraph(strcat(path,'Paper/Figures/Robustness_beforeGM'), 'Quarter',  '', 1.4);

%%
%------ AFTER the great moderation : 1985-2014 ------%
%------ Sample Selection ------%
% SVAR
DataSVAR = [oil gdp emp pcepi ff];
DataSVAR = DataSVAR(105:end,:);
nt = size(DataSVAR,1);         
ns = size(DataSVAR,2);           
bw_bw = 100; 
for is = 1:ns; 
    tmp = bw_trend(DataSVAR(:,is),bw_bw);
   	DataSVAR_trend_robust2b(:,is)= tmp;
   	DataSVAR(:,is) = DataSVAR(:,is) - DataSVAR_trend_robust2b(:,is); 
end;
xmean = nanmean(DataSVAR)';                                       
mult = sqrt((sum(~isnan(DataSVAR))-1)./sum(~isnan(DataSVAR)));     
xstd = (nanstd(DataSVAR).*mult)';                                  
DATA_SVAR_robust2b = (DataSVAR - repmat(xmean',nt,1))./repmat(xstd',nt,1);  
% FAVAR
DATA_FACTORS = data_factors(105:end,:);
F = factorsFAVAR(DATA_FACTORS,numFt);
DATA_FAVAR_robust2b = [DATA_SVAR_robust2b F];
%------ Estimation ------%
% SVAR
%p = hqcL;                  
constant = 1;             
hor = 20;                
iter = 500;                
conf = [90 60];            
Cf_SVAR_robust2b = ir(DATA_SVAR_robust2b,p,constant,hor,'cholimpact');
irfs_SVAR = permute(Cf_SVAR_robust1,[3 1 2]);
irfs_SVAR = irfs_SVAR(:,:,1);
% FAVAR
Cf_FAVAR_robust2b = ir(DATA_FAVAR_robust2b,p,constant,hor,'cholimpact');
[ULbf_FAVAR_robust2b,~,~,~,ULb2f_FAVAR_robust2b] = bootbands(DATA_FAVAR_robust2b,p,constant,iter,hor,conf,'cholimpact');
irfs_FAVAR = permute(Cf_FAVAR_robust2b,[3 1 2]);
irfs_FAVAR = irfs_FAVAR(:,:,1);
conf_FAVAR1 = permute(ULbf_FAVAR_robust2b,[3 1 2 4]);
conf_FAVAR_down1 = conf_FAVAR1(:,:,1,1);
conf_FAVAR_up1 = conf_FAVAR1(:,:,1,2);
conf_FAVAR2 = permute(ULb2f_FAVAR_robust2b,[3 1 2 4]);
conf_FAVAR_down2 = conf_FAVAR2(:,:,1,1);
conf_FAVAR_up2 = conf_FAVAR2(:,:,1,2);
plot_titles = { 'Oil Price', 'Real GDP' , 'Total Employment', 'Inflation'};
figure(10);
for ii = 1:4
subplot(2,2,ii); 
plot(irfs_SVAR(:,ii),'red','LineWidth',1.5,'LineStyle','-');hold on;     
plot(irfs_FAVAR(:,ii),'black','LineWidth',2);title(plot_titles(:,ii)); 
formatGraph('Quarter', '');
hold on;plot(xlim, [0 0], 'color', 0.5*[1 1 1],'LineWidth',1.1);                     
hold on;plot(conf_FAVAR_down1(:,ii),'black','LineStyle',':');hold on;plot(conf_FAVAR_up1(:,ii),'black','LineStyle',':')
hold on;plot(conf_FAVAR_down2(:,ii),'black','LineStyle','-');hold on;plot(conf_FAVAR_up2(:,ii),'black','LineStyle','-')
end
saveGraph(strcat(path,'Paper/Figures/Robustness_afterGM'), 'Quarter',  '', 1.4);


%%
%------------------------------------------------------------------------%
%------ Lags in the VAR : 3 - 8 ------%
%------------------------------------------------------------------------%

%------ 3 lags ------%
% SVAR
p = 3;                  
constant = 1;             
hor = 20;                
iter = 500;                
conf = [90 60];            
Cf_SVAR_3lags = ir(DATA_SVAR,p,constant,hor,'cholimpact');
irfs_SVAR = permute(Cf_SVAR_3lags,[3 1 2]);
irfs_SVAR = irfs_SVAR(:,:,1);
% FAVAR
Cf_FAVAR_3lags = ir(DATA_FAVAR,p,constant,hor,'cholimpact');
[ULbf_FAVAR_3lags,~,~,~,ULb2f_FAVAR_3lags] = bootbands(DATA_FAVAR,p,constant,iter,hor,conf,'cholimpact');
irfs_FAVAR = permute(Cf_FAVAR_3lags,[3 1 2]);
irfs_FAVAR = irfs_FAVAR(:,:,1);
conf_FAVAR1 = permute(ULbf_FAVAR_3lags,[3 1 2 4]);
conf_FAVAR_down1 = conf_FAVAR1(:,:,1,1);
conf_FAVAR_up1 = conf_FAVAR1(:,:,1,2);
conf_FAVAR2 = permute(ULb2f_FAVAR_3lags,[3 1 2 4]);
conf_FAVAR_down2 = conf_FAVAR2(:,:,1,1);
conf_FAVAR_up2 = conf_FAVAR2(:,:,1,2);
plot_titles = { 'Oil Price', 'Real GDP' , 'Total Employment', 'Inflation'};
figure(11);
for ii = 1:4
subplot(2,2,ii); 
plot(irfs_SVAR(:,ii),'red','LineWidth',1.5,'LineStyle','-');hold on;      
plot(irfs_FAVAR(:,ii),'black','LineWidth',2);title(plot_titles(:,ii));  
formatGraph('Quarter', '');
hold on;plot(xlim, [0 0], 'color', 0.5*[1 1 1],'LineWidth',1.1);                     
hold on;plot(conf_FAVAR_down1(:,ii),'black','LineStyle',':');hold on;plot(conf_FAVAR_up1(:,ii),'black','LineStyle',':')
hold on;plot(conf_FAVAR_down2(:,ii),'black','LineStyle','-');hold on;plot(conf_FAVAR_up2(:,ii),'black','LineStyle','-')
end
saveGraph(strcat(path,'Paper/Figures/Robustness_3lags'), 'Quarter',  '', 1.4);

%------ SVAR 2/3 lags ------%
% SVAR
p = 3;                  
constant = 1;             
hor = 20;                
iter = 500;                
conf = [90 60];            
Cf_SVAR_23lags = ir(DATA_SVAR,p,constant,hor,'cholimpact');
irfs_SVAR = permute(Cf_SVAR_23lags,[3 1 2]);
irfs_SVAR = irfs_SVAR(:,:,1);
% FAVAR
Cf_FAVAR_3lags = ir(DATA_FAVAR,p,constant,hor,'cholimpact');
[ULbf_FAVAR_3lags,~,~,~,ULb2f_FAVAR_3lags] = bootbands(DATA_FAVAR,p,constant,iter,hor,conf,'cholimpact');
irfs_FAVAR = permute(Cf_FAVAR_3lags,[3 1 2]);
irfs_FAVAR = irfs_FAVAR(:,:,1);
conf_FAVAR1 = permute(ULbf_FAVAR_3lags,[3 1 2 4]);
conf_FAVAR_down1 = conf_FAVAR1(:,:,1,1);
conf_FAVAR_up1 = conf_FAVAR1(:,:,1,2);
conf_FAVAR2 = permute(ULb2f_FAVAR_3lags,[3 1 2 4]);
conf_FAVAR_down2 = conf_FAVAR2(:,:,1,1);
conf_FAVAR_up2 = conf_FAVAR2(:,:,1,2);
plot_titles = { 'Oil Price', 'Real GDP' , 'Total Employment', 'Inflation'};
figure(11);
for ii = 1:4
subplot(2,2,ii); 
plot(irfs_SVAR(:,ii),'red','LineWidth',1.5,'LineStyle','-');hold on;      
plot(irfs_FAVAR(:,ii),'black','LineWidth',2);title(plot_titles(:,ii));  
formatGraph('Quarter', '');
hold on;plot(xlim, [0 0], 'color', 0.5*[1 1 1],'LineWidth',1.1);                     
hold on;plot(conf_FAVAR_down1(:,ii),'black','LineStyle',':');hold on;plot(conf_FAVAR_up1(:,ii),'black','LineStyle',':')
hold on;plot(conf_FAVAR_down2(:,ii),'black','LineStyle','-');hold on;plot(conf_FAVAR_up2(:,ii),'black','LineStyle','-')
end
saveGraph(strcat(path,'Paper/Figures/Robustness_3lags'), 'Quarter',  '', 1.4);

%%
%------ 8 lags ------%
% SVAR
p = 8;                  
constant = 1;             
hor = 20;                
iter = 500;                
conf = [90 60];            
Cf_SVAR_8lags = ir(DATA_SVAR,p,constant,hor,'cholimpact');
irfs_SVAR = permute(Cf_SVAR_8lags,[3 1 2]);
irfs_SVAR = irfs_SVAR(:,:,1);
% FAVAR
Cf_FAVAR_8lags = ir(DATA_FAVAR,p,constant,hor,'cholimpact');
[ULbf_FAVAR_8lags,~,~,~,ULb2f_FAVAR_8lags] = bootbands(DATA_FAVAR,p,constant,iter,hor,conf,'cholimpact');
irfs_FAVAR = permute(Cf_FAVAR_8lags,[3 1 2]);
irfs_FAVAR = irfs_FAVAR(:,:,1);
conf_FAVAR1 = permute(ULbf_FAVAR_8lags,[3 1 2 4]);
conf_FAVAR_down1 = conf_FAVAR1(:,:,1,1);
conf_FAVAR_up1 = conf_FAVAR1(:,:,1,2);
conf_FAVAR2 = permute(ULb2f_FAVAR_8lags,[3 1 2 4]);
conf_FAVAR_down2 = conf_FAVAR2(:,:,1,1);
conf_FAVAR_up2 = conf_FAVAR2(:,:,1,2);
plot_titles = { 'Oil Price', 'Real GDP' , 'Total Employment', 'Inflation'};
figure(12);
for ii = 1:4
subplot(2,2,ii); 
plot(irfs_SVAR(:,ii),'red','LineWidth',1.5,'LineStyle','-');hold on;      
plot(irfs_FAVAR(:,ii),'black','LineWidth',2);title(plot_titles(:,ii));  
formatGraph('Quarter', '');
hold on;plot(xlim, [0 0], 'color', 0.5*[1 1 1],'LineWidth',1.1);                     
hold on;plot(conf_FAVAR_down1(:,ii),'black','LineStyle',':');hold on;plot(conf_FAVAR_up1(:,ii),'black','LineStyle',':')
hold on;plot(conf_FAVAR_down2(:,ii),'black','LineStyle','-');hold on;plot(conf_FAVAR_up2(:,ii),'black','LineStyle','-')
end
saveGraph(strcat(path,'Paper/Figures/Robustness_8lags'), 'Quarter',  '', 1.4);


%%
%------------------------------------------------------------------------%
%------ Factors in the FAVAR : 3 ------%
%------------------------------------------------------------------------%

%------ 3 factors ------%
% SVAR
p = 2;                  
constant = 1;             
hor = 20;                
iter = 500;                
conf = [90 60];            
Cf_SVAR_3factors = ir(DATA_SVAR,p,constant,hor,'cholimpact');
irfs_SVAR = permute(Cf_SVAR_3factors,[3 1 2]);
irfs_SVAR = irfs_SVAR(:,:,1);
% FAVAR
numFt = 8;
F = factorsFAVAR(data_factors,numFt); % set the appropriate number of factors
DATA_FAVAR = [DATA_SVAR F];
Cf_FAVAR_3factors = ir(DATA_FAVAR,p,constant,hor,'cholimpact');
[ULbf_FAVAR_3factors,~,~,~,ULb2f_FAVAR_3factors] = bootbands(DATA_FAVAR,p,constant,iter,hor,conf,'cholimpact');
irfs_FAVAR = permute(Cf_FAVAR_3factors,[3 1 2]);
irfs_FAVAR = irfs_FAVAR(:,:,1);
conf_FAVAR1 = permute(ULbf_FAVAR_3factors,[3 1 2 4]);
conf_FAVAR_down1 = conf_FAVAR1(:,:,1,1);
conf_FAVAR_up1 = conf_FAVAR1(:,:,1,2);
conf_FAVAR2 = permute(ULb2f_FAVAR_3factors,[3 1 2 4]);
conf_FAVAR_down2 = conf_FAVAR2(:,:,1,1);
conf_FAVAR_up2 = conf_FAVAR2(:,:,1,2);
plot_titles = { 'Oil Price', 'Real GDP' , 'Total Employment', 'Inflation'};
figure(13);
for ii = 1:4
subplot(2,2,ii); 
plot(irfs_SVAR(:,ii),'red','LineWidth',1.5,'LineStyle','-');hold on;      
plot(irfs_FAVAR(:,ii),'black','LineWidth',2);title(plot_titles(:,ii));  
formatGraph('Quarter', '');
hold on;plot(xlim, [0 0], 'color', 0.5*[1 1 1],'LineWidth',1.1);                     
hold on;plot(conf_FAVAR_down1(:,ii),'black','LineStyle',':');hold on;plot(conf_FAVAR_up1(:,ii),'black','LineStyle',':')
hold on;plot(conf_FAVAR_down2(:,ii),'black','LineStyle','-');hold on;plot(conf_FAVAR_up2(:,ii),'black','LineStyle','-')
end
saveGraph(strcat(path,'Paper/Figures/Robustness_3factors'), 'Quarter',  '', 1.4);

%%
%------------------------------------------------------------------------%
%------ FAVAR vs SVAR vs lags ------%
%------------------------------------------------------------------------%

for p=2:3
    %------ SVAR estimation ------%

    % Parameters        
    constant = 1;              % include the constant
    hor = 20;                  % horizon
    iter = 500;                % number of iterations
    conf = [90 60];            % level for confidence bands

    % Estimation of the SVAR
    Cf_SVAR = ir(DATA_SVAR,p,constant,hor,'cholimpact');
    [ULbf_SVAR,~,~,~,ULb2f_SVAR] = bootbands(DATA_SVAR,p,constant,iter,hor,conf,'cholimpact');

    % Plot resutls
    irfs_SVAR = permute(Cf_SVAR,[3 1 2]);
    irfs_SVAR = irfs_SVAR(:,:,1);
    conf_SVAR1 = permute(ULbf_SVAR,[3 1 2 4]);
    conf_SVAR_down1 = conf_SVAR1(:,:,1,1);
    conf_SVAR_up1 = conf_SVAR1(:,:,1,2);
    conf_SVAR2 = permute(ULb2f_SVAR,[3 1 2 4]);
    conf_SVAR_down2 = conf_SVAR2(:,:,1,1);
    conf_SVAR_up2 = conf_SVAR2(:,:,1,2);
    plot_titles = { 'Oil Price', 'Real GDP' , 'Total Employment', 'Inflation'};

    %------ FAVAR estimation ------%

    F = factorsFAVAR(data_factors, 2);

    % Data
    DATA_FAVAR = [DATA_SVAR F];

    % Parameters
    %p = CVrel;                  % number of lags
    constant = 1;              % include the constant
    hor = 20;                  % horizon
    iter = 500;                % number of iterations
    conf = [90 60];            % level for confidence bands

    % Estimation of the FAVAR
    Cf_FAVAR = ir(DATA_FAVAR,p,constant,hor,'cholimpact');
    [ULbf_FAVAR,~,~,~,ULb2f_FAVAR] = bootbands(DATA_FAVAR,p,constant,iter,hor,conf,'cholimpact');

    % Plot resutls
    irfs_FAVAR = permute(Cf_FAVAR,[3 1 2]);
    irfs_FAVAR = irfs_FAVAR(:,:,1);
    conf_FAVAR1 = permute(ULbf_FAVAR,[3 1 2 4]);
    conf_FAVAR_down1 = conf_FAVAR1(:,:,1,1);
    conf_FAVAR_up1 = conf_FAVAR1(:,:,1,2);
    conf_FAVAR2 = permute(ULb2f_FAVAR,[3 1 2 4]);
    conf_FAVAR_down2 = conf_FAVAR2(:,:,1,1);
    conf_FAVAR_up2 = conf_FAVAR2(:,:,1,2);
    plot_titles = { 'Oil Price', 'Real GDP' , 'Total Employment', 'Inflation'};

    for ii = 1:4
    subplot(2,2,ii); 
    % IRFs of the SVAR
    plot(irfs_SVAR(:,ii),'black','LineWidth',1.5,'LineStyle','-');hold on; 
    % IRFs of the FAVAR
    plot(irfs_FAVAR(:,ii),'red','LineWidth',2);title(plot_titles(:,ii));  
    formatGraph('Quarter', '');
    % Horizontal line at 0
    hold on;plot(xlim, [0 0], 'color', 0.5*[1 1 1],'LineWidth',1.1);                     
    % 90% confidence interval
    hold on;plot(conf_SVAR_down1(:,ii),'black','LineStyle',':');hold on;plot(conf_SVAR_up1(:,ii),'black','LineStyle',':')
    % 60% confidence interval
    hold on;plot(conf_SVAR_down2(:,ii),'black','LineStyle','-');hold on;plot(conf_SVAR_up2(:,ii),'black','LineStyle','-')
    end
    legend('SVAR with errors','FAVAR', 'Location','southeast')
    saveGraph(strcat(path,strcat('Paper/Figures/rob_SVAR_', num2str(p))), 'Quarter',  '', 1.4);
    close all;
    for ii = 1:4
    subplot(2,2,ii); 
    % IRFs of the SVAR
    plot(irfs_SVAR(:,ii),'red','LineWidth',1.5,'LineStyle','-');hold on; 
    % IRFs of the FAVAR
    plot(irfs_FAVAR(:,ii),'black','LineWidth',2);title(plot_titles(:,ii));  
    formatGraph('Quarter', '');
    % Horizontal line at 0
    hold on;plot(xlim, [0 0], 'color', 0.5*[1 1 1],'LineWidth',1.1);                     
    % 90% confidence interval
    hold on;plot(conf_FAVAR_down1(:,ii),'black','LineStyle',':');hold on;plot(conf_FAVAR_up1(:,ii),'black','LineStyle',':')
    % 60% confidence interval
    hold on;plot(conf_FAVAR_down2(:,ii),'black','LineStyle','-');hold on;plot(conf_FAVAR_up2(:,ii),'black','LineStyle','-')
    end
    legend('SVAR','FAVAR with errors', 'Location','southeast')
    saveGraph(strcat(path,strcat('Paper/Figures/rob_FAVAR_', num2str(p))), 'Quarter',  '', 1.4);
    close all;
end
