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
[fit_Htop stats]=fit(x,y,'power1');
IV_tph_BE=(IV_tph_ERA+IV_tph_NOAA)/2-IV_ventalt;
ytph=y./IV_tph_BE(mask);
masktph=ytph<=1.25;
[fit_Htoptph statstph]=fit(x(masktph),y(masktph),'power1')

xlist=logspace(min(log10(x)),max(log10(x)),1000)
%https://uk.mathworks.com/help/curvefit/predint.html
fit_int= predint(fit_Htop,xlist,0.95,'functional');
fit_int2= predint(fit_Htop,xlist,0.95,'observation');
% fit_int= predint(fit_Htop,xlist);
figure(2)
subplot(1,2,1)

X=[xlist,fliplr(xlist)];                %#create continuous x value array for plotting
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];              %#create y values for out and then back
huncert=fill(X,Y,'r','FaceColor','#003f5c','EdgeColor','none','FaceAlpha',0.25); 
hold on
huncert2=plot(xlist,fit_int2(:,1),'--','Color','#003f5c','LineWidth',1)
plot(xlist,fit_int2(:,2),'--','Color','#003f5c','LineWidth',1)
hfit=plot(xlist,fit_Htop(xlist),'-','Color','#003f5c','LineWidth',3)


hdata=errorbar(x,y,dy,dy,dx_l,dx_u,'ok','MarkerFaceColor','k','MarkerSize',3,'CapSize',1,'LineWidth',0.03)

set(gca,'XScale','log')



ci_Htop = confint(fit_Htop);
ci_Htop=0.5*(ci_Htop(2,:)-ci_Htop(1,:));
c_Htop = coeffvalues(fit_Htop);
ci_Htop=round(ci_Htop,3);c_Htop=round(c_Htop,3);
cvaltext=strcat('All events: R^2=',num2str(round(stats.rsquare,2)));
cvaltexttph=strcat('Tropospheric/lower stratospheric: R^2=',num2str(round(statstph.rsquare,2)));


% legend('IVESPA','Best fit','Mastin et al. (2009)','Sparks et al. (1997)','Wilson and Walker (1987)')

xlabel('$\rm \overline{MER} \ (kg \ s^{-1})$','Interpreter','Latex')
ylabel('$\rm \overline{H}_{top}$ (km a.v.l.)','Interpreter','Latex')
xlim([0.5*min(x) 2*max(x)])
ylim([0 42])
text(20,0.7*42,cvaltext,'FontSize',10)
text(20,0.66*42,cvaltexttph,'FontSize',10)
text(-0.08,0.98,'a)','Units','normalized','FontWeight','Bold','FontSize',12)



%==========================================================================
%Height vs MER (top height)
%==========================================================================
%DRE density = 2500
IV_tph_BE=(IV_tph_ERA+IV_tph_NOAA)/2-IV_ventalt;
x=IV_TEM_BE./IV_duration_BE;
dx_l=x.*((IV_TEM_UL./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
dx_u=x.*((IV_TEM_UU./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
y=IV_Htop_BE./IV_tph_BE;
dy=IV_Htop_U./IV_tph_BE;
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dx_l=dx_l(mask);dy=dy(mask);dx_u=dx_u(mask);
[fit_Htop stats]=fit(x,y,'power1')
masktph=y<=1.25;
[fit_Htoptph statstph]=fit(x(masktph),y(masktph),'power1')
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


hdata=errorbar(x,y,dy,dy,dx_l,dx_u,'ok','MarkerFaceColor','k','MarkerSize',3,'CapSize',1,'LineWidth',0.03)

set(gca,'XScale','log')



ci_Htop = confint(fit_Htop);
ci_Htop=0.5*(ci_Htop(2,:)-ci_Htop(1,:));
c_Htop = coeffvalues(fit_Htop);
ci_Htop=round(ci_Htop,3);c_Htop=round(c_Htop,3);
cvaltext=strcat('All events: R^2=',num2str(round(stats.rsquare,2)));
cvaltexttph=strcat('Tropospheric/lower stratospheric: R^2=',num2str(round(statstph.rsquare,2)));
% legend([hdata hfit huncert],'IVESPA',cvaltext,'95% confidence interval on fit function','Location','Northwest')


% legend('IVESPA','Best fit','Mastin et al. (2009)','Sparks et al. (1997)','Wilson and Walker (1987)')

xlabel('$\rm \overline{MER} \ (kg \ s^{-1})$','Interpreter','Latex')
ylabel('$\rm \overline{H}_{top}$/tropopause height','Interpreter','Latex')
xlim([0.5*min(x) 2*max(x)])
ylim([0 8])
text(20,0.7*8,cvaltext,'FontSize',10)
text(20,0.66*8,cvaltexttph,'FontSize',10)
text(-0.08,0.98,'b)','Units','normalized','FontWeight','Bold','FontSize',12)