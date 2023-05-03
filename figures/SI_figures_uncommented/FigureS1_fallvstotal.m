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

%==========================================================================
%Height vs total MER (top height)
%==========================================================================
%DRE density = 2500
%again I use generic variables x and y so that I can easily copy this bloc
%of code and change variable definition in only one place

%I replace unknown PDC masses by 0 (assuming there was no PDC; debatable)
IV_PDC_BE(isnan(IV_PDC_BE))=0;

%x is the total mass eruption rate, i.e. the sum of the fall+PDC mass
%divided by duration
x=(IV_TEM_BE+IV_PDC_BE)./IV_duration_BE;

%I try to calculate MER uncertainty using regular propagation rules, but
%with non-symmetric uncertainty for mass that's likely shady... I end up
%applying the regular propagation rules but with distinct upper and lower
%bounds
dx_l=x.*(((IV_TEM_UL.^2+IV_PDC_U.^2).^0.5./(IV_TEM_BE+IV_PDC_BE)).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
dx_u=x.*(((IV_TEM_UU.^2+IV_PDC_U.^2).^0.5./(IV_TEM_BE+IV_PDC_BE)).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;

%y is the top height
y=IV_Htop_BE;
dy=IV_Htop_U;

%here I remove NaN values
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dx_l=dx_l(mask);dy=dy(mask);dx_u=dx_u(mask);

%Here I fit a power law and then calculate the uncertainty related to the
%fitting procedure
[fit_Htoptot stats]=fit(x,y,'power1')
xlist=logspace(min(log10(x)),max(log10(x)),1000)
%https://uk.mathworks.com/help/curvefit/predint.html
fit_int= predint(fit_Htoptot,xlist,0.95,'functional');
fit_int2= predint(fit_Htoptot,xlist,0.95,'observation');
% fit_int= predint(fit_Htop,xlist);

%I create a figure with six panels and plot in the first one
figure(2)
subplot(1,2,1)


%plotting the data and the fit with its uncertainty...
X=[xlist,fliplr(xlist)];                %#create continuous x value array for plotting
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];              %#create y values for out and then back
huncert=fill(X,Y,'r','FaceColor','#003f5c','EdgeColor','none','FaceAlpha',0.25); 
hold on
hfit=plot(xlist,fit_Htoptot(xlist),'-','Color','#003f5c','LineWidth',3)
huncert2=plot(xlist,fit_int2(:,1),'--','Color','#003f5c','LineWidth',1)
plot(xlist,fit_int2(:,2),'--','Color','#003f5c','LineWidth',1)

maskpdc=IV_PDC_BE(mask)./(IV_TEM_BE(mask)+IV_PDC_BE(mask))<0.1;
hdata=errorbar(x(maskpdc),y(maskpdc),dy(maskpdc),dy(maskpdc),dx_l(maskpdc),dx_u(maskpdc),'^k','MarkerFaceColor',[1 1 1],'MarkerSize',5,'CapSize',1,'LineWidth',0.03)
maskpdc=IV_PDC_BE(mask)./(IV_TEM_BE(mask)+IV_PDC_BE(mask))>=0.1 & IV_PDC_BE(mask)./(IV_TEM_BE(mask)+IV_PDC_BE(mask))<0.5 ;
hdata2=errorbar(x(maskpdc),y(maskpdc),dy(maskpdc),dy(maskpdc),dx_l(maskpdc),dx_u(maskpdc),'sk','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',5,'CapSize',1,'LineWidth',0.03)
maskpdc=IV_PDC_BE(mask)./(IV_TEM_BE(mask)+IV_PDC_BE(mask))>=0.5 ;
hdata3=errorbar(x(maskpdc),y(maskpdc),dy(maskpdc),dy(maskpdc),dx_l(maskpdc),dx_u(maskpdc),'ko','MarkerFaceColor','k','MarkerSize',5,'CapSize',1,'LineWidth',0.03)


set(gca,'XScale','log')



%retrieving fit parameters including confidence intervals

ci_Htoptot = confint(fit_Htoptot);
ci_Htoptot=0.5*(ci_Htoptot(2,:)-ci_Htoptot(1,:));
c_Htoptot = coeffvalues(fit_Htoptot);
ci_Htoptot=round(ci_Htoptot,3);c_Htoptot=round(c_Htoptot,3);
cvaltext={strcat('Fit parameters: a=',num2str(c_Htoptot(1)),'\pm',num2str(ci_Htoptot(1)),';');strcat('b=',num2str(c_Htoptot(2)),'\pm',num2str(ci_Htoptot(2)),'; R^2=',num2str(round(stats.rsquare,2)))}

%legend, label, etc
% legend('IVESPA','Best fit','Mastin et al. (2009)','Sparks et al. (1997)','Wilson and Walker (1987)')
xlabel('$\rm \overline{MER} \ (kg \ s^{-1}, total)$','Interpreter','Latex')

ylabel('$\rm \overline{H}_{top}$ (km a.v.l.)','Interpreter','Latex')
xlim([0.5*min(x) 2*max(x)])
ylim([0 42])

text(100,25,cvaltext)
legend([hdata hdata2 hdata3 hfit huncert huncert2],'$\rm PDC \ mass \ \% \ < \ 0.1$','$\rm 0.1 \ \leq \ PDC \ mass \ \% \ < \ 0.5$','$\rm 0.5 \ \leq \ PDC \ mass \ \% $','Power law fit ($\rm \overline{H}_{top} = a \times \overline{MER}^b$)','$\rm Functionnal \ uncertainty \ (95\% \ confidence \ level)$','$\rm Prediction \ uncertainty \ (95\% \ confidence \ level)$','Location','Northwest','Interpreter','Latex')

title('a) $\rm \overline{MER}$ based on fallout+PDC mass','Interpreter','Latex')
%Then I'm repeating the same thing 4 times for other heights or for the
%fall mass eruption rate (no pdc included). The last panels plots all the
%fit together, without IVESPA data, for comparison
%==========================================================================
%Height vs MER (top height)
%==========================================================================
%DRE density = 2500
x=IV_TEM_BE./IV_duration_BE;
dx_l=x.*((IV_TEM_UL./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
dx_u=x.*((IV_TEM_UU./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
y=IV_Htop_BE;
dy=IV_Htop_U;
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dx_l=dx_l(mask);dy=dy(mask);dx_u=dx_u(mask);
[fit_Htop stats]=fit(x,y,'power1')
xlist=logspace(min(log10(x)),max(log10(x)),1000)
%https://uk.mathworks.com/help/curvefit/predint.html
fit_int= predint(fit_Htop,xlist,0.95,'functional');
fit_int2= predint(fit_Htop,xlist,0.95,'observation');
% fit_int= predint(fit_Htop,xlist);
figure(2)
subplot(1,2,2)

X=[xlist,fliplr(xlist)];                %#create continuous x value array for plotting
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];              %#create y values for out and then back
huncert=fill(X,Y,'r','FaceColor','#003f5c','EdgeColor','none','FaceAlpha',0.25); 
hold on
huncert2=plot(xlist,fit_int2(:,1),'--','Color','#003f5c','LineWidth',1)
plot(xlist,fit_int2(:,2),'--','Color','#003f5c','LineWidth',1)
hfit=plot(xlist,fit_Htop(xlist),'-','Color','#003f5c','LineWidth',3)

maskpdc=IV_PDC_BE(mask)./(IV_TEM_BE(mask)+IV_PDC_BE(mask))<0.1;
hdata=errorbar(x(maskpdc),y(maskpdc),dy(maskpdc),dy(maskpdc),dx_l(maskpdc),dx_u(maskpdc),'^k','MarkerFaceColor',[1 1 1],'MarkerSize',5,'CapSize',1,'LineWidth',0.03)
maskpdc=IV_PDC_BE(mask)./(IV_TEM_BE(mask)+IV_PDC_BE(mask))>=0.1 & IV_PDC_BE(mask)./(IV_TEM_BE(mask)+IV_PDC_BE(mask))<0.5 ;
hdata2=errorbar(x(maskpdc),y(maskpdc),dy(maskpdc),dy(maskpdc),dx_l(maskpdc),dx_u(maskpdc),'sk','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',5,'CapSize',1,'LineWidth',0.03)
maskpdc=IV_PDC_BE(mask)./(IV_TEM_BE(mask)+IV_PDC_BE(mask))>=0.5 ;
hdata3=errorbar(x(maskpdc),y(maskpdc),dy(maskpdc),dy(maskpdc),dx_l(maskpdc),dx_u(maskpdc),'ko','MarkerFaceColor','k','MarkerSize',5,'CapSize',1,'LineWidth',0.03)

set(gca,'XScale','log')



ci_Htop = confint(fit_Htop);
ci_Htop=0.5*(ci_Htop(2,:)-ci_Htop(1,:));
c_Htop = coeffvalues(fit_Htop);
ci_Htop=round(ci_Htop,3);c_Htop=round(c_Htop,3);
cvaltext={strcat('Fit parameters: a=',num2str(c_Htop(1)),'\pm',num2str(ci_Htop(1)),';');strcat('b=',num2str(c_Htop(2)),'\pm',num2str(ci_Htop(2)),'; R^2=',num2str(round(stats.rsquare,2)))}
% legend([hdata hfit huncert],'IVESPA',cvaltext,'95% confidence interval on fit function','Location','Northwest')


% legend('IVESPA','Best fit','Mastin et al. (2009)','Sparks et al. (1997)','Wilson and Walker (1987)')

xlabel('$\rm \overline{MER} \ (kg \ s^{-1}, fallout)$','Interpreter','Latex')
ylabel('$\rm \overline{H}_{top}$ (km a.v.l.)','Interpreter','Latex')
xlim([0.5*min(x) 2*max(x)])
ylim([0 42])
text(100,25,cvaltext)
title('b) $\rm \overline{MER}$ based on fallout mass only','Interpreter','Latex')