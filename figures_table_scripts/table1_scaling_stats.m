% 	Written by Thomas J. Aubry, May 2023.
% 	Department of Earth and Environmental Sciences, University of Exeter
%   E-mail: t.aubry@exeter.ac.uk
% 	Please cite the corresponding paper if you use this script

clear
close all

addpath('../functions/') %to access useful functions
addpath('../') %to access data
load_IVESPA %scripts reading IVESPA data


IV_MER_BE=IV_TEM_BE./IV_duration_BE;%define MER
y=IV_Htop_BE;%define y as top height
mask=~isnan(IV_MER_BE) & ~isnan(IV_Htop_BE);%create a mask for events with values for top height

%=======================================================================
%Create different sets of weight
%=======================================================================
%1-uniform weight (i.e. no weight)
w_none=ones(size(IV_Htop_BE(mask)));

%2-give the same weight to each eruption
eru=IV_eruption(mask);
w_eruption=NaN(size(eru));
for iu=1:length(eru)
    w_eruption(iu)=1/sum(eru==eru(iu));%i.e. weight inversely proportional to number of event in each eruption
end


%3-weight according to uncertainty on predicted minus observed height

%first, for events with no TEM uncertainty, assume a factor of three
%uncertainty to get lower and upper bount uncertainty
dTEM_UL=IV_TEM_UL;dTEM_UL(isnan(dTEM_UL))=2*IV_TEM_BE(isnan(dTEM_UL))/3;
dTEM_UU=IV_TEM_UU;dTEM_UU(isnan(dTEM_UU))=2*IV_TEM_BE(isnan(dTEM_UU));
%calculate TEM relative uncertainty as the mean relative lower and upper
%bound uncertainty
dTEM=(dTEM_UU./IV_TEM_BE+dTEM_UL./IV_TEM_BE)/2;
%use standard propagation rules to get a MER uncertainty
IV_MER_U=IV_MER_BE.*(dTEM.^2+(IV_duration_U./IV_duration_BE).^2).^0.5;

%calculate uncertainty on height predicted using the power law scaling
%using standard propagation rules.
IV_Hpred=0.3448*IV_MER_BE.^0.226;
IV_Hpred_U=IV_Hpred.*(IV_MER_U./IV_MER_BE)*0.226;

%Define weight as inversely proportional to the uncertainty on the
%difference between predicted and observed height
w_uncertainty=IV_Hpred_U(mask).^2+IV_Htop_U(mask).^2;
w_uncertainty=1./w_uncertainty;


%4-weight according to value of interpretation flags for MER (duration and
%TEM) and top height
w_flag=1+IV_Htop_BE_flag+max(IV_TEM_BE_flag,IV_duration_BE_flag);
w_flag=1./w_flag(mask);

%5-weight according to the previous three factors considered (i.e. simply
%the product of these three weights)
w_all=w_flag.*w_uncertainty.*w_eruption;

%Create a matrix containing all 5 sets of weights. Note that Matlab will
%automatically standardize weights so that their sum is equal to 1.
weights=[w_none w_eruption w_uncertainty w_flag w_all];


%=======================================================================
%Calibrate model and calculate performance metrics
%=======================================================================
adjR2=NaN(8,5);%preallocate space for adjusted R2
AIC=NaN(8,5);%preallocate space for bias-corrected AIC
coef_a=NaN(7,5);%preallocate space for various scaling empirical parameters
coef_b=NaN(7,5);
coef_c=NaN(7,5);
coef_d=NaN(7,5);

%for each set of weight, we will repeat the calibration/evaluation process
%for all 8 scalings considered
for iw=1:5
    
w=weights(:,iw); %select correct set of weight



%model 1 = empirical power law
im=2;
x=[IV_MER_BE];%define matrix containing independent variables (only MER in this case)
start=[1 1];%define starting condition
modelFun = @(b,x) scaling1(b,x);%define function containing the scaling equation (see scaling1 function file)
wnlm = fitnlm(x(mask,:),y(mask),modelFun,start,'Weight',w);%fit scaling using Matlab's non-linear model fitting function
adjR2(im,iw)=wnlm.Rsquared.Adjusted;%save adjusted R2
%save fit parameters (2 in this case, prefactor a and exponent b)
coef_a(im,iw)=wnlm.Coefficients.Estimate(1);
coef_b(im,iw)=wnlm.Coefficients.Estimate(2);
AIC(im,iw)=wnlm.ModelCriterion.AICc;%save bias-corrected AIC


%model 2 = empirical power law with wind and stratification
im=3;
x=[IV_MER_BE IV_N_BE IV_W_BE];%this time independent variables include MER, wind and stratification
start=[1 1 1 1];
modelFun = @(b,x) scaling2(b,x);%see scaling2 function for details
wnlm = fitnlm(x(mask,:),y(mask),modelFun,start,'Weight',w);
adjR2(im,iw)=wnlm.Rsquared.Adjusted;
%save the four fit parameters
coef_a(im,iw)=wnlm.Coefficients.Estimate(1);
coef_b(im,iw)=wnlm.Coefficients.Estimate(2);
coef_c(im,iw)=wnlm.Coefficients.Estimate(3);
coef_d(im,iw)=wnlm.Coefficients.Estimate(4);
AIC(im,iw)=wnlm.ModelCriterion.AICc;

%model 3 =MTT
im=4;
x=[IV_MER_BE IV_N_BE];%this time independent variables include MER and stratification
start=[1];
modelFun = @(b,x) scaling3(b,x);%see scaling3 function for details
wnlm = fitnlm(x(mask,:),y(mask),modelFun,start,'Weight',w);
adjR2(im,iw)=wnlm.Rsquared.Adjusted;
coef_a(im,iw)=wnlm.Coefficients.Estimate(1);
AIC(im,iw)=wnlm.ModelCriterion.AICc;

%model 4 =Hewett et al 71
im=5;
x=[IV_MER_BE IV_N_BE IV_W_BE];%independent variables: MER, wind and stratification

start=[1];
modelFun = @(b,x) scaling4(b,x);%see scaling4 function for details
wnlm = fitnlm(x(mask,:),y(mask),modelFun,start,'Weight',w);
adjR2(im,iw)=wnlm.Rsquared.Adjusted;
coef_a(im,iw)=wnlm.Coefficients.Estimate(1);
AIC(im,iw)=wnlm.ModelCriterion.AICc;

%model 5 = Degruyter and Bonadonna 2012; b is beta^2/alpha^1.5
im=6;
x=[IV_MER_BE IV_N_BE IV_Vs_BE];%independent variables: MER, stratification, and wind regime parameter Vs
start=[1 1];
modelFun = @(b,x) scaling5(b,x);%see scaling5 function for details
wnlm = fitnlm(x(mask,:),y(mask),modelFun,start,'Weight',w);
adjR2(im,iw)=wnlm.Rsquared.Adjusted;
coef_a(im,iw)=wnlm.Coefficients.Estimate(1);
coef_b(im,iw)=wnlm.Coefficients.Estimate(2);
AIC(im,iw)=wnlm.ModelCriterion.AICc;


%model 6 = Woodhouse et al. 2013; b is beta/alpha
im=7;
x=[IV_MER_BE IV_N_BE IV_Riw_BE];%independent variables: MER, stratification, and wind regime parameter Riw
start=[1 1];
modelFun = @(b,x) scaling6(b,x);%see scaling6 function for details
wnlm = fitnlm(x(mask,:),y(mask),modelFun,start,'Weight',w);
adjR2(im,iw)=wnlm.Rsquared.Adjusted;
coef_a(im,iw)=wnlm.Coefficients.Estimate(1);
coef_b(im,iw)=wnlm.Coefficients.Estimate(2);
AIC(im,iw)=wnlm.ModelCriterion.AICc;


%model 7 = Aubry al. 2017; b is beta/alpha

im=8;
x=[IV_MER_BE IV_N_BE IV_Ws_BE];%independent variables: MER, stratification, and wind regime parameter WS
start=[1 1];
modelFun = @(b,x) scaling7(b,x);%see scaling7 function for details
wnlm = fitnlm(x(mask,:),y(mask),modelFun,start,'Weight',w);
adjR2(im,iw)=wnlm.Rsquared.Adjusted;
coef_a(im,iw)=wnlm.Coefficients.Estimate(1);
coef_b(im,iw)=wnlm.Coefficients.Estimate(2);
AIC(im,iw)=wnlm.ModelCriterion.AICc;


% % %model 1 = Mastin et al. (2009)
im=1;
x=[IV_MER_BE];
start=[1 1];
%In this specific case, we do not calibrate the model but instead directly
%use the coefficients from Mastin et al. (2009)
modelFun = @(x) scaling1([0.304 0.241],x);
modelpred=modelFun(x);

%the adjusted R2 is calculated "manually" using the sum of squares of errors (SSE) and
%total sum of square (SST)
SSE=sum((-(modelpred(mask)-y(mask)).*(w.^0.5)).^2);
SST=sum((-(sum(y(mask).*w)/sum(w)-y(mask)).*(w.^0.5)).^2);
rsqr=1-SSE/SST;
adjR2(im,iw)=1-(1-rsqr)*(length(y(mask))-1)/(length(y(mask))-size(x,2)-1);


end


%=======================================================================
%create a table with parameter values
%=======================================================================

%the code below just concatenate together all the R2/parameter values
%calculated above to facilitate integration in a text table.
tabval=strings(8,5);

for iw=1:5
    for im=1:8
        
        txt=strcat('R^2=',num2str(adjR2(im,iw),2));
        if im>1
           txt=strcat(txt,', a=',num2str(coef_a(im,iw),2));  
        end
if im~=4 & im~=5 & im~=1
    txt=strcat(txt,', b=',num2str(coef_b(im,iw),2));
end

if im==3
    txt=strcat(txt,', c=',num2str(coef_c(im,iw),2),', d=',num2str(coef_d(im,iw),2));
end
tabval(im,iw)=txt;
    end
end

