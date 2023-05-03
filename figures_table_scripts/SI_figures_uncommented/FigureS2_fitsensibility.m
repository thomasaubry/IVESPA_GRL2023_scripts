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

x=IV_TEM_BE./IV_duration_BE;
dx_l=x.*((IV_TEM_UL./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
dx_u=x.*((IV_TEM_UU./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
y=IV_Htop_BE;
dy=IV_Htop_U;
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dx_l=dx_l(mask);dy=dy(mask);dx_u=dx_u(mask);
IV_year=IV_year(mask);
IV_MER_BE=IV_MER_BE(mask);
reg=IV_region(mask);
IV_ivid=IV_ivid(mask);
reglist=unique(reg);
IV_volcano=IV_volcano(mask);
IV_W_BE=IV_W_BE(mask);
IV_N_BE=IV_N_BE(mask);
IV_heighttech=IV_heighttech(mask);
IV_style=IV_style(mask);
IV_PI_BE=IV_PI_BE(mask);
%==========================================================================
%Height vs MER (top height)
%==========================================================================
%DRE density = 2500

[fit_Htop stats]=fit(x,y,'power1')
fit_Htop_ref=fit_Htop;
xlist=logspace(min(log10(x)),max(log10(x)),1000)
%https://uk.mathworks.com/help/curvefit/predint.html
fit_int= predint(fit_Htop,xlist,0.95,'functional');

% fit_int= predint(fit_Htop,xlist);
figure(2)




%==========================================================================
%Height vs MER (top height)
%==========================================================================
%DRE density = 2500
% mask=IV_year<1990;
% [fit_Htop stats]=fit(x(mask),y(mask),'power1')
% 
% figure(2)
% subplot(1,3,1)
% hfit=plot(xlist,fit_Htop_ref(xlist),':','Color','k','LineWidth',2)
% hold on
% hfit=plot(xlist,fit_Htop(xlist),'-','Color','#003f5c','LineWidth',2)
% hold on
% hdata=errorbar(x(mask),y(mask),dy(mask),dy(mask),dx_l(mask),dx_u(mask),'ko','MarkerFaceColor','#003f5c','MarkerSize',3,'CapSize',1,'LineWidth',0.03)
% 
% 
% mask=IV_year>=1990 & IV_year<2006;
% [fit_Htop stats]=fit(x(mask),y(mask),'power1')
% hfit=plot(xlist,fit_Htop(xlist),'-','Color','#bc5090','LineWidth',2)
% hold on
% hdata=errorbar(x(mask),y(mask),dy(mask),dy(mask),dx_l(mask),dx_u(mask),'ks','MarkerFaceColor','#bc5090','MarkerSize',3,'CapSize',1,'LineWidth',0.03)
% 
% 
% 
% mask=IV_year>=2006;
% [fit_Htop stats]=fit(x(mask),y(mask),'power1')
% hfit=plot(xlist,fit_Htop(xlist),'-','Color','#ffa600','LineWidth',2)
% hold on
% hdata=errorbar(x(mask),y(mask),dy(mask),dy(mask),dx_l(mask),dx_u(mask),'k^','MarkerFaceColor','#ffa600','MarkerSize',3,'CapSize',1,'LineWidth',0.03)
% 
% legend('Fit for all events','Fit for 1900-1989 events','1900-1989 events','Fit for 1990-2005 events','1990-2005 events','Fit for 2006-2016 events','2006-2016 events')
% 
% 
% set(gca,'XScale','log')
% 
% xlabel('$\rm \overline{MER} \ (kg \ s^{-1})$','Interpreter','Latex')
% ylabel('$\rm \overline{H}_{top}$ (km a.v.l.)','Interpreter','Latex')
% xlim([0.5*min(x) 2*max(x)])
% ylim([0 max(y)*1.1])
% 
% title('a) Sensitivity to event year')

%==========================================================================
%Height vs MER (top height)
%==========================================================================
%DRE density = 2500
mask=IV_MER_BE>quantile(IV_MER_BE,0.1);
[fit_Htop stats]=fit(x(mask),y(mask),'power1')
mask=IV_MER_BE<=quantile(IV_MER_BE,0.1);
figure(2)
subplot(1,2,1)
hfit=plot(xlist,fit_Htop_ref(xlist),':','Color','k','LineWidth',2)
hold on
hfit=plot(xlist,fit_Htop(xlist),'-','Color','#003f5c','LineWidth',2)
hold on
hdata=errorbar(x(mask),y(mask),dy(mask),dy(mask),dx_l(mask),dx_u(mask),'ko','MarkerFaceColor','#003f5c','MarkerSize',3,'CapSize',1,'LineWidth',0.03)


mask=IV_MER_BE<quantile(IV_MER_BE,0.9);
[fit_Htop stats]=fit(x(mask),y(mask),'power1')
mask=IV_MER_BE>=quantile(IV_MER_BE,0.1);
hfit=plot(xlist,fit_Htop(xlist),'-','Color','#bc5090','LineWidth',2)
hold on
hdata=errorbar(x(mask),y(mask),dy(mask),dy(mask),dx_l(mask),dx_u(mask),'ks','MarkerFaceColor','#bc5090','MarkerSize',3,'CapSize',1,'LineWidth',0.03)



mask=IV_MER_BE<=quantile(IV_MER_BE,0.9) & IV_MER_BE>=quantile(IV_MER_BE,0.1);
[fit_Htop stats]=fit(x(mask),y(mask),'power1')
hfit=plot(xlist,fit_Htop(xlist),'-','Color','#ffa600','LineWidth',2)
hold on

hdata=errorbar(x(mask),y(mask),dy(mask),dy(mask),dx_l(mask),dx_u(mask),'k^','MarkerFaceColor','#ffa600','MarkerSize',3,'CapSize',1,'LineWidth',0.03)

legend('Fit for all events','Fit without 10\% of events with lowest $\rm \overline{MER}$','10\% of events with lowest $\rm \overline{MER}$','Fit without 10\% of events with highest $\rm \overline{MER}$','10\% of events with highest $\rm \overline{MER}$','Fit for events within \newline the $\rm 10^{th}$ and $\rm 90^{th}$ $\rm \overline{MER}$ quantile','Events within the $\rm 10^{th}$ and $\rm 90^{th}$ $\rm \overline{MER}$ percentile','Interpreter','Latex')


set(gca,'XScale','log')

xlabel('$\rm \overline{MER} \ (kg \ s^{-1})$','Interpreter','Latex')
ylabel('$\rm \overline{H}_{top}$ (km a.v.l.)','Interpreter','Latex')
xlim([0.5*min(x) 2*max(x)])
ylim([0 max(y)*1.1])

title('a) Sensitivity to $\rm \overline{MER}$ range','Interpreter','Latex')






%==========================================================================
%Height vs MER (top height)
%==========================================================================
%DRE density = 2500
mask=~strcmp(IV_ivid,'STM1902_01') & ~strcmp(IV_ivid,'QUI1932_01') & ~strcmp(IV_ivid,'BEZ1956_01');
[fit_Htop stats]=fit(x(mask),y(mask),'power1')
figure(2)
subplot(1,2,2)
hfit=plot(xlist,fit_Htop_ref(xlist),':','Color','k','LineWidth',2)
hold on
hfit=plot(xlist,fit_Htop(xlist),'-','Color','#003f5c','LineWidth',2)
hold on
hdata=errorbar(x(~mask),y(~mask),dy(~mask),dy(~mask),dx_l(~mask),dx_u(~mask),'ko','MarkerFaceColor','#003f5c','MarkerSize',3,'CapSize',1,'LineWidth',0.03)


mask=~strcmp(IV_ivid,'COT2015_04');
[fit_Htop stats]=fit(x(mask),y(mask),'power1')

hfit=plot(xlist,fit_Htop(xlist),'-','Color','#7a5195','LineWidth',2)
hold on
hdata=errorbar(x(~mask),y(~mask),dy(~mask),dy(~mask),dx_l(~mask),dx_u(~mask),'ks','MarkerFaceColor','#7a5195','MarkerSize',3,'CapSize',1,'LineWidth',0.03)



mask=IV_volcano~=313030;
[fit_Htop stats]=fit(x(mask),y(mask),'power1')
hfit=plot(xlist,fit_Htop(xlist),'-','Color','#ffa600','LineWidth',2)
hold on

hdata=errorbar(x(~mask),y(~mask),dy(~mask),dy(~mask),dx_l(~mask),dx_u(~mask),'k^','MarkerFaceColor','#ffa600','MarkerSize',3,'CapSize',1,'LineWidth',0.03)


mask=IV_volcano~=211060;
[fit_Htop stats]=fit(x(mask),y(mask),'power1')
hfit=plot(xlist,fit_Htop(xlist),'-','Color','#ef5675','LineWidth',2)
hold on

hdata=errorbar(x(~mask),y(~mask),dy(~mask),dy(~mask),dx_l(~mask),dx_u(~mask),'kp','MarkerFaceColor','#ef5675','MarkerSize',3,'CapSize',1,'LineWidth',0.03)



legend('Fit for all events','Fit without SM1902, Quiz1932, Bez1956','SM1902, Quiz1932, Bez1956','Fit without Cotopaxi 2015 phase 4','Cotopaxi 2015 phase 4','Fit without Redoubt events','Redoubt events','Fit without Etna events','Etna events')

set(gca,'XScale','log')

xlabel('$\rm \overline{MER} \ (kg \ s^{-1})$','Interpreter','Latex')
ylabel('$\rm \overline{H}_{top}$ (km a.v.l.)','Interpreter','Latex')
xlim([0.5*min(x) 2*max(x)])
ylim([0 max(y)*1.1])

title('$\rm b) \ Sensitivity \ to \ inclusion \ of \ specific \ events$','Interpreter','Latex')









