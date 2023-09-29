function STATS=kwtest(x)
%Kruskal-Wallis test for the non parametric ANOVA
%In statistics, the Kruskal–Wallis one-way analysis of variance by ranks (named
%after William Kruskal and W. Allen Wallis) is a non-parametric method for
%testing equality of population medians among groups. It is identical to a
%one-way analysis of variance with the data replaced by their ranks. It is an
%extension of the Mann–Whitney U test to 3 or more groups. Since it is a
%non-parametric method, the Kruskal–Wallis test does not assume a normal
%population, unlike the analogous one-way analysis of variance. However, the
%test does assume an identically-shaped and scaled distribution for each group,
%except for any difference in medians. The exact distribution of the
%Kruskal-Wallis statistics is very time, space and memory expensive, so it is
%approximated using other distribution.
%The MatLab function KRUSKALWALLIS only uses the chi-square distribution that
%is the most conservative (and this means that it accepts the H0 hypothesis more
%than you want). This function computes also the F, Beta and Gamma
%approximations (the F distribution is the less conservative and this means that
%it refuses the H0 more than you want), giving a more informative output (if you
%want, see http://www.jmu.edu/assessment/JPM%20AERA%20SP%2008.pdf). If you
%believe that the p-value you choose is smaller than your cut-off (usually
%0.05), you can use my function Dunn-Sidak to isolate the differences among
%groups (http://www.mathworks.com/matlabcentral/fileexchange/12827).
%
% Syntax: 	STATS=kwtest(X)
%      
%     Inputs:
%           X - data matrix (Size of matrix must be n-by-2; data=column 1, group=column 2). 
%
%     Outputs:
%           - Statistics of each group (samples, median, sum of ranks and mean rank)
%           - Correction factor for ties and H statistics
%           - Chi-square approximation
%           - F approximation
%           - Beta approximation
%           - Gamma approximation
%        If STATS nargout was specified the results will be stored in the STATS
%        struct.
%
%      Example: 
%
%                             Sample
%                   -------------------------
%                       1       2        3    
%                   -------------------------
%                      7.89    8.84     8.65
%                      9.16    9.92    10.70
%                      7.34    7.20    10.24
%                     10.28    9.25     8.62
%                      9.12    9.45     9.94
%                      9.24    9.14    10.55
%                      8.40    9.99    10.13
%                      8.60    9.21     9.78          
%                      8.04    9.06     9.01
%                      8.45
%                      9.51
%                      8.15
%                      7.69
%                   -------------------------
%
%       Data matrix must be:
% x=[7.79 9.16 7.64 10.28 9.12 9.24 8.40 8.60 8.04 8.45 9.51 8.15 7.69 ...
% 8.84 9.92 7.20 9.25 9.45 9.14 9.99 9.21 9.06 8.65 10.70 10.24 8.62 ...
% 9.94 10.55 10.13 9.78 9.01; ...
% repmat(1,1,13) repmat(2,1,9) repmat(3,1,9)]'
%
%           Calling on Matlab the function: kwtest(x)
%
%           Answer is:
%
% KRUSKAL-WALLIS TEST
% --------------------------------------------------------------------------------
%     Group    Samples    Median    Ranks_sum    Mean_rank
%     _____    _______    ______    _________    _________
% 
%     1        13         8.45      146          11.231   
%     2         9         9.21      152          16.889   
%     3         9         9.94      198              22   
% 
% Correction factor for ties: 0.0000	H: 7.5823
% --------------------------------------------------------------------------------
%  
% Chi-square approximation (the most conservative)
% --------------------------------------------------------------------------------
%     Chi_square    df    two_tailed_p_value
%     __________    __    __________________
% 
%     7.5823        2     0.02257           
% 
% F-statistic approximation (the less conservative)
% --------------------------------------------------------------------------------
%       F       df_num    df_denom    two_tailed_p_value
%     ______    ______    ________    __________________
% 
%     5.0734    2         29          0.013473          
% 
% Beta distribution approximation
% --------------------------------------------------------------------------------
%        B       alpha      beta     two_tailed_p_value
%     _______    ______    ______    __________________
% 
%     0.28779    0.9438    11.489    0.018072          
% 
% Gamma distribution approximation
% --------------------------------------------------------------------------------
%       G       alpha      beta     two_tailed_p_value
%     ______    ______    ______    __________________
% 
%     7.5823    1.1035    1.8124    0.019             
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2009). KWTEST: Kruskal-Wallis non parametric test for ANalysis Of VAriance
% http://www.mathworks.com/matlabcentral/fileexchange/25860

%Input Error handling
p = inputParser;
%addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','nonempty','ncols',2}));
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'real','nonempty','ncols',2}));
parse(p,x);
assert(all(x(:,2) == fix(x(:,2))),'Warning: all elements of column 2 of input matrix must be whole numbers')
clear p


tr=repmat('-',1,80);% set divisor
disp('KRUSKAL-WALLIS TEST')
disp(tr)
N=size(x,1);%total observation
%One first ranks all the observations without regard for which treatment group
%they are in.
[R,T]=tiedrank(x(:,1));

%Next, compute the ranks sum for each group.
k=max(x(:,2)); %number of groups
KW=zeros(4,k); %matrix preallocation
for I=1:k
   KW(1,I)=length(x(x(:,2)==I)); %elements of the i-th group
   KW(2,I)=median(x(x(:,2)==I),1); %median of the i-th group
   KW(3,I)=sum(R(x(:,2)==I)); %the ranks sum
   KW(4,I)=KW(3,I)/KW(1,I); %the mean rank
end
for I=1:k
   KW(1,I)=length(rmmissing(x(x(:,2)==I))); %elements of the i-th group
   KW(2,I)=nanmedian(x(x(:,2)==I),1); %median of the i-th group
   KW(3,I)=nansum(R(x(:,2)==I)); %the ranks sum
   KW(4,I)=KW(3,I)/KW(1,I); %the mean rank
end
disp(table((1:k)',KW(1,:)',KW(2,:)',KW(3,:)',KW(4,:)','VariableNames',{'Group','Samples','Median','Ranks_sum','Mean_rank'}))
%the mean rank computed when tratment has no effect.
Rbar=(N+1)/2;
%We will use the sum of squared deviations between each sample group's average
%rank and the overall average rank, weighted by the sizes of each group, as a
%measure of variability between the observations and what you would expect if
%the hypothesis of no treatment effect was true. Call this sum D.  
D=sum(KW(1,:).*(KW(4,:)-Rbar).^2);
Hbiased=12*D/N/(N+1); %Kruskal-Wallis statistic uncorrected for ties
if T==0 %if there are not ties
    CF=0;
    H=Hbiased;    
else
    CF=1-2*T/N/(N^2-1); %correction factor for ties
    H=Hbiased/CF;
end
fprintf('Correction factor for ties: %0.4f\tH: %0.4f\n',CF,H)
disp(tr); disp(' ')

df=k-1; %degrees of freedom
fprintf('Chi-square approximation (the most conservative)\n')
disp(tr)
P1=1-chi2cdf(H,df);  %probability associated to the Chi-squared-statistic.
disp(table(H,df,P1,'VariableNames',{'Chi_square','df','two_tailed_p_value'}))

fprintf('F-statistic approximation (the less conservative)\n')
disp(tr)
dfd=N-df; %degrees of freedom for denominator
F=((dfd+1)*H)/((k-1)*(N-1-H));  %F-statistic approximation.
P2=1-fcdf(F,df,N-k-1);  %probability associated to the F-statistic.
disp(table(F,df,dfd,P2,'VariableNames',{'F','df_num','df_denom','two_tailed_p_value'}))

m=k-1; %the expected value
s2=2*df-2*(3*k^2-6*k+N*(2*k^2-6*k+1))/(5*N*(N+1))-(6/5*sum(1./KW(1,:))); %variance
eta=(N^3-sum(KW(1,:).^3))/N/(N+1); %maximum value of H

fprintf('Beta distribution approximation\n')
disp(tr)
B=H/eta;
alpha1=m*((m*(eta-m)-s2)/(eta*s2));
beta1=alpha1*((eta-m)/m);
P3=1-betacdf(B,alpha1,beta1);  %probability associated to the Beta distribution.
disp(table(B,alpha1,beta1,P3,'VariableNames',{'B','alpha','beta','two_tailed_p_value'}))

fprintf('Gamma distribution approximation\n')
disp(tr)
alpha2=m^2/s2;
beta2=s2/m;
P4=1-gamcdf(H,alpha2,beta2);  %probability associated to the Gamma Distribution.
disp(table(H,alpha2,beta2,P4,'VariableNames',{'G','alpha','beta','two_tailed_p_value'}))

if nargout
    STATS.Hbiased=Hbiased;
    STATS.CF=CF;
    STATS.H=H;
    STATS.Chi_square.chi2=H;
    STATS.Chi_square.df=df;
    STATS.Chi_square.pvalue=P1;
    STATS.F.F=F;
    STATS.F.dfn=df;
    STATS.F.dfd=dfd;
    STATS.F.pvalue=P2;
    STATS.Beta.m=m;
    STATS.Beta.s2=s2;
    STATS.Beta.eta=eta;
    STATS.Beta.B=B;
    STATS.Beta.alpha=alpha1;
    STATS.Beta.beta=beta1;
    STATS.Beta.pvalue=P3;
    STATS.Gamma.m=m;
    STATS.Gamma.s2=s2;
    STATS.Gamma.G=H;
    STATS.Gamma.alpha=alpha2;
    STATS.Gamma.beta=beta2;
    STATS.Gamma.pvalue=P4;
end
