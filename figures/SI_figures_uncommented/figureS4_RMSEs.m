% 	Written by Thomas J. Aubry, May 2023.
% 	Department of Earth and Environmental Sciences, University of Exeter
%   E-mail: t.aubry@exeter.ac.uk
% 	Please cite the corresponding paper if you use this script
%   Apologies for the lack of comments in SI figure script! Feel free to
%   email me for help using this script.


clear
close all

addpath('../functions/')
addpath('../')
load_IVESPA

x=IV_MER_BE;
dx_l=x.*((IV_TEM_UL./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
dx_u=x.*((IV_TEM_UU./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
y=IV_Htop_BE;
dy=IV_Htop_U;
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dx_l=dx_l(mask);dy=dy(mask);dx_u=dx_u(mask);
IV_year=IV_year(mask);
reg=IV_region(mask);
IV_ivid=IV_ivid(mask);
reglist=unique(reg);
IV_volcano=IV_volcano(mask);
IV_W_BE=IV_W_BE(mask);
IV_N_BE=IV_N_BE(mask);
IV_heighttech=IV_heighttech(mask);
IV_style=IV_style(mask);
IV_PI_BE=IV_PI_BE(mask);
IV_duration_BE=IV_duration_BE(mask);
IV_TEM_BE=IV_TEM_BE(mask);

fontlab=5;
%==========================================================================
%Height predictions from MER
%==========================================================================
%DRE density = 2500
load('fit_Mtop.mat')
load('fit_Htop.mat')
c_Htop = coeffvalues(fit_Htop);
c_Mtop = coeffvalues(fit_Mtop);

ylist=linspace(min(y),max(y),1000);


figure(2)
subplot(2,2,1)


hfit_H=plot(y,fit_Htop(x),'o','Color','#003f5c','MarkerFaceColor','#003f5c','MarkerSize',3)
RMSE_H=round(sqrt(mean(((fit_Htop(x)-y)).^2)),1);
RMSRE_H=round(100*sqrt(mean(((fit_Htop(x)-y)./y).^2)),1);
hold on
hfit_M=plot(y,10.^((log10(x)-c_Mtop(2))/c_Mtop(1)),'^','Color','#003f5c','MarkerFaceColor',hex2rgb('#003f5c')+0.75*([1 1 1]-hex2rgb('#003f5c')),'MarkerSize',3)
RMSE_M=round(sqrt(mean(((10.^((log10(x)-c_Mtop(2))/c_Mtop(1))-y)).^2)),1);
RMSRE_M=round(100*sqrt(mean(((10.^((log10(x)-c_Mtop(2))/c_Mtop(1))-y)./y).^2)),1);

hfit_mast=plot(y,2*(x/2500).^0.241,'ks','MarkerFaceColor','r','MarkerSize',3)
RMSE_mast=round(sqrt(mean(((2*(x/2500).^0.241-y)).^2)),1);
RMSRE_mast=round(100*sqrt(mean(((2*(x/2500).^0.241-y)./y).^2)),1);

hline=plot([0 38],[0 38],'k:','LineWidth',1.5)
xlim([0 38])
ylim([0 38])
lg_H=['$\rm \overline{H}_{top}=0.345 \times \overline{MER}^{0.226},$' newline strcat('$\rm RMSE=',num2str(RMSE_H),'\ km \ (',num2str(RMSRE_H),'\%$)')];
lg_M=['$\rm log(\overline{MER})=2.83+3.54 \times log(\overline{H}_{top}),$' newline strcat('$\rm RMSE=',num2str(RMSE_M),'\ km \ (',num2str(RMSRE_M),'\%$)')];
lg_mast=['Mastin et al. (2009),' newline strcat('$\rm RMSE=',num2str(RMSE_mast),'\ km \ (',num2str(RMSRE_mast),'\%$)')];
legend(lg_H,lg_M,lg_mast,'Interpreter','Latex')
legend('boxoff')

xlabel('Observed $\rm \overline{H}_{top}$ (km a.v.l.)','Interpreter','Latex')
ylabel('Predicted $\rm \overline{H}_{top}$ (km a.v.l.)','Interpreter','Latex')
title('a) $\rm \overline{H}_{top} \ from \ \overline{MER}$','Interpreter','Latex')
set(gca,'Xtick',0:5:35,'Ytick',0:5:35)

% errmetr=((fit_Htop(x)-y)./y).^2;
% [val ind]=sort(errmetr);
% for i=0:4
%   text(y(ind(end-i))+0.75,fit_Htop(x(ind(end-i)))+0.75,IV_ivid(ind(end-i)),'Interpreter','none','FontSize',fontlab);   
% end
% errmetr=((fit_Htop(x)-y)).^2;
% [val ind]=sort(errmetr);
% for i=0:4
%   text(y(ind(end-i))+0.75,fit_Htop(x(ind(end-i)))+0.75,IV_ivid(ind(end-i)),'Interpreter','none','FontSize',fontlab);   
% end


subplot(2,2,2)

hfit_H=plot(x,(y/c_Htop(1)).^(1/c_Htop(2)),'o','Color','#003f5c','MarkerFaceColor','#003f5c','MarkerSize',3)
RMSE_H=round(sqrt(mean(((log10((y/c_Htop(1)).^(1/c_Htop(2)))-log10(x))).^2)),2);

hold on
hfit_M=plot(x,10.^fit_Mtop(log10(y)),'^','Color','#003f5c','MarkerFaceColor',hex2rgb('#003f5c')+0.75*([1 1 1]-hex2rgb('#003f5c')),'MarkerSize',3)
RMSE_M=round(sqrt(mean((fit_Mtop(log10(y))-log10(x)).^2)),2);

set(gca,'Xscale','log','Yscale','log')
hfit_mast=plot(x,2500*(y/2).^(1/0.241),'ks','MarkerFaceColor','r','MarkerSize',3)
RMSE_mast=round(sqrt(mean((log10(2500*(y/2).^(1/0.241))-log10(x)).^2)),2);


hline=plot([10 2*10^9],[10 2*10^9],'k:','LineWidth',1.5)
xlim([10 2*10^9])
ylim([10 2*10^9])
lg_H=strcat('$\rm \overline{H}_{top}=0.345\times \overline{MER}^{0.226}$, lRMSE=',num2str(RMSE_H));
lg_M=strcat('$\rm log(\overline{MER})=2.83+3.54\times log(\overline{H}_{top})$, lRMSE=',num2str(RMSE_M));
lg_mast=strcat('Mastin et al. (2009), lRMSE=',num2str(RMSE_mast));
legend(lg_H,lg_M,lg_mast,'Interpreter','Latex')
legend('boxoff')

xlabel('Observed $\rm \overline{MER}\ (kg\ s^{-1})$','Interpreter','Latex')
ylabel('Predicted $\rm \overline{MER}\ (kg\ s^{-1})$','Interpreter','Latex')
title('b) $\rm \overline{MER} \ from \ \overline{H}_{top}$','Interpreter','Latex')
set(gca,'Xtick',10.^(1:9),'Ytick',10.^(1:9))


% errmetr=(fit_Mtop(log10(y))-log10(x)).^2;
% [val ind]=sort(errmetr);
% for i=0:4
%   text(x(ind(end-i))*1.4,10.^fit_Mtop(log10(y(ind(end-i))))*1.4,IV_ivid(ind(end-i)),'Interpreter','none','FontSize',fontlab);   
% end



subplot(2,2,3)

hfit_H=plot(IV_duration_BE,IV_TEM_BE./((y/c_Htop(1)).^(1/c_Htop(2))),'o','Color','#003f5c','MarkerFaceColor','#003f5c','MarkerSize',3)
RMSE_H=round(sqrt(mean(((log10(IV_TEM_BE./(((y/c_Htop(1)).^(1/c_Htop(2)))))-log10(IV_duration_BE))).^2)),2);

hold on
hfit_M=plot(IV_duration_BE,IV_TEM_BE./(10.^fit_Mtop(log10(y))),'^','Color','#003f5c','MarkerFaceColor',hex2rgb('#003f5c')+0.75*([1 1 1]-hex2rgb('#003f5c')),'MarkerSize',3)
RMSE_M=round(sqrt(mean((log10(IV_TEM_BE./10.^fit_Mtop(log10(y)))-log10(IV_duration_BE)).^2)),2);

set(gca,'Xscale','log','Yscale','log')
hfit_mast=plot(IV_duration_BE,IV_TEM_BE./(2500*(y/2).^(1/0.241)),'ks','MarkerFaceColor','r','MarkerSize',3)
RMSE_mast=round(sqrt(mean((log10(IV_TEM_BE./(2500*(y/2).^(1/0.241)))-log10(IV_duration_BE)).^2)),2);


hline=plot([10 3*10^7],[10 3*10^7],'k:','LineWidth',1.5)
xlim([10 3*10^7])
ylim([10 3*10^7])
lg_H=strcat('$\rm \overline{H}_{top}=0.345 \times \overline{MER}^{0.226}$, lRMSE=',num2str(RMSE_H));
lg_M=strcat('$\rm log(\overline{MER})=2.83+3.54 \times log(\overline{H}_{top})$, lRMSE=',num2str(RMSE_M));
lg_mast=strcat('Mastin et al. (2009): lRMSE=',num2str(RMSE_mast));
legend(lg_H,lg_M,lg_mast,'Interpreter','Latex')
legend('boxoff')

xlabel('$\rm Observed \ duration \ (s)$','Interpreter','Latex')
ylabel('$\rm Predicted \ duration \ (s)$','Interpreter','Latex')
title('c) Duration from TEM and $\rm \overline{H}_{top}$','Interpreter','Latex')
set(gca,'Xtick',10.^(1:7),'Ytick',10.^(1:7))


% errmetr=(log10(IV_TEM_BE./10.^fit_Mtop(log10(y)))-log10(IV_duration_BE)).^2;
% [val ind]=sort(errmetr);
% for i=0:4
%   text(IV_duration_BE(ind(end-i))*1.4,IV_TEM_BE(ind(end-i))./(10.^fit_Mtop(log10(y(ind(end-i)))))*1.4,IV_ivid(ind(end-i)),'Interpreter','none','FontSize',fontlab);   
% end



subplot(2,2,4)

hfit_H=plot(IV_TEM_BE,(y/c_Htop(1)).^(1/c_Htop(2)).*IV_duration_BE,'o','Color','#003f5c','MarkerFaceColor','#003f5c','MarkerSize',3)
RMSE_H=round(sqrt(mean(((log10((y/c_Htop(1)).^(1/c_Htop(2)).*IV_duration_BE)-log10(IV_TEM_BE))).^2)),2);

hold on
hfit_M=plot(IV_TEM_BE,10.^fit_Mtop(log10(y)).*IV_duration_BE,'^','Color','#003f5c','MarkerFaceColor',hex2rgb('#003f5c')+0.75*([1 1 1]-hex2rgb('#003f5c')),'MarkerSize',3)
RMSE_M=round(sqrt(mean((log10(10.^fit_Mtop(log10(y)).*IV_duration_BE)-log10(IV_TEM_BE)).^2)),2);

set(gca,'Xscale','log','Yscale','log')
hfit_mast=plot(IV_TEM_BE,2500*(y/2).^(1/0.241).*IV_duration_BE,'ks','MarkerFaceColor','r','MarkerSize',3)
RMSE_mast=round(sqrt(mean((log10(2500*(y/2).^(1/0.241).*IV_duration_BE)-log10(IV_TEM_BE)).^2)),2);


hline=plot([2*10^5 2*10^13],[2*10^5 2*10^13],'k:','LineWidth',1.5)
xlim([2*10^5 2*10^13])
ylim([2*10^5 2*10^13])
lg_H=strcat('$\rm \overline{H}_{top}=0.345 \times \overline{MER}^{0.226}$, lRMSE=',num2str(RMSE_H));
lg_M=strcat('$\rm log(\overline{MER})=2.83+3.54 \times log(\overline{H}_{top})$, lRMSE=',num2str(RMSE_M));
lg_mast=strcat('Mastin et al. (2009): lRMSE=',num2str(RMSE_mast));
legend(lg_H,lg_M,lg_mast,'Interpreter','Latex')
legend('boxoff')

xlabel('$\rm Observed \ tephra \ fallout \ mass \ (kg)$','Interpreter','Latex')
ylabel('$\rm Predicted \ tephra \ fallout \ mass \ (kg)$','Interpreter','Latex')
title('d) TEM from $\rm \overline{H}_{top}$ and duration','Interpreter','Latex')
set(gca,'Xtick',10.^(5:13),'Ytick',10.^(5:13))

% errmetr=(log10(10.^fit_Mtop(log10(y)).*IV_duration_BE)-log10(IV_TEM_BE)).^2;
% [val ind]=sort(errmetr);
% for i=0:4
%   text(IV_TEM_BE(ind(end-i))*1.4,10.^fit_Mtop(log10(y(ind(end-i)))).*IV_duration_BE(ind(end-i))*1.4,IV_ivid(ind(end-i)),'Interpreter','none','FontSize',fontlab);   
% end


