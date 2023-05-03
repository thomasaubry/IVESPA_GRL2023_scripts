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

eventlist={'RED1990_06','STR2003_01','RED1990_07','REV2002_01','NEV1985_01','RED1990_02','QUI1932_01','STM1902_01','ETN2013_01','ANA2003_02','EYJ2010_04','CHA2008_03','SHM1997_01','MER2010_01','CHA2008_01','BEZ1956_01'};

eventlist1={'MER2010_01','CHA2008_01','BEZ1956_01','MIY2000_01','COT2015_01'};
eventlist2={'MER2010_01','CHA2008_01','MIY2000_01','COT2015_01','RED1990_05'};
eventlist3={'MER2010_01','CHA2008_01','MIY2000_01','COT2015_01','GRI2004_02','NEV1985_01'};
eventlist4={'MER2010_01','CHA2008_01','MIY2000_01','COT2015_01','ETN2016_02'};


%==========================================================================
%Height vs MER (top height)
%==========================================================================
%DRE density = 2500
x=IV_MER_BE;
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
% fit_int= predint(fit_Htop,xlist);
z=IV_morpho(mask);
mask_weak=strcmp('weak',z);
mask_strong=strcmp('strong',z);
mask_unknown=~mask_weak & ~mask_strong;
meanwind=IV_W_BE(mask);
IV_heighttech=IV_heighttech(mask);


% IV_PI_BE=IV_PI_BE(mask);

fontlab=5;
meroffset=1.2;meroffset2=1.1
hoffset=1;
%==========================================================================
%res plot 1
%==========================================================================
hres=IV_Htop_BE./fit_Htop(IV_MER_BE);
x=IV_W_BE./IV_W_ref;

figure
subplot(2,2,2)
mask=abs(IV_latitude)<25;
plot(x(mask),hres(mask),'ks','MarkerSize',5,'MarkerFaceColor',[0.6 0 0])
hold on
mask=abs(IV_latitude)>25;
plot(x(mask),hres(mask),'k>','MarkerSize',5,'MarkerFaceColor',[0.6 0.8 1])
hold on
%ylim([0.2 3.2])
%xlim([2 39])
ylabel('$\rm \overline{H}_{top}^{std}$','interpreter','latex')  
xlabel('$\rm Detrended\ wind \ speed$','Interpreter','Latex')


mask=~isnan(hres);
[fitwind stats]=fit(x(mask),hres(mask),'poly1')
[cc pp]=corrcoef(x(mask),hres(mask));cc=round(cc(1,2),2);pp(1,2)
plot(linspace(min(x),max(x),500),fitwind(linspace(min(x),max(x),500)),':','LineWidth',2,'Color',[0 0 0])
legend('Tropical events','Extra-tropical events',strcat('Linear trend: \bf r=',num2str(cc)),'Location','NorthEast')
xlim([0.9*min(x) 1.1*max(x)])
ylim([min(hres)-0.15 max(hres)+0.15])




% for jj=1:length(eventlist2)
% ii=strcmp(IV_ivid,eventlist2{jj});text(x(ii)+1,hres(ii),eventlist2{jj},'Interpreter','none','FontSize',fontlab);
% end


%==========================================================================
%res plot 1
%==========================================================================
hres=IV_Htop_BE./fit_Htop(IV_MER_BE);
%hres=hres./fit_wshear(IV_Wshear_BE);
x=IV_N_BE./IV_N_ref;


subplot(2,2,1)

mask=abs(IV_latitude)<25;
plot(x(mask),hres(mask),'ks','MarkerSize',5,'MarkerFaceColor',[0.6 0 0])
hold on
mask=abs(IV_latitude)>25;
plot(x(mask),hres(mask),'k>','MarkerSize',5,'MarkerFaceColor',[0.6 0.8 1])

hold on
% ylim([-8 18])
% xlim([2 39])
ylabel('$\rm \overline{H}_{top}^{std}$','interpreter','latex')  
xlabel('$\rm Detrended\ Brunt$ V\"ais\"al\"a frequency','Interpreter','Latex')


mask=~isnan(hres);
[fitwind stats]=fit(x(mask),hres(mask),'poly1')
[cc pp]=corrcoef(x(mask),hres(mask));cc=round(cc(1,2),2);pp(1,2)
plot(linspace(min(x),max(x),500),fitwind(linspace(min(x),max(x),500)),':','LineWidth',2,'Color',[0 0 0])
% legend('Magmatic event','Phreatomagmatic','Phreatic','Unknown Style','Linear fit (R^2=0)','Location','NorthEast')
legend('Tropical events','Extra-tropical events',strcat('Linear trend: r=',num2str(cc)),'Location','NorthEast')
xlim([0.95*min(x) max(x)+0.05*min(x)])
ylim([min(hres)-0.15 max(hres)+0.15])



% for jj=1:length(eventlist1)
% ii=strcmp(IV_ivid,eventlist1{jj});text(x(ii)+0.0003,hres(ii),eventlist1{jj},'Interpreter','none','FontSize',fontlab);
% end
%set(gca,'Yscale','log','Xscale','log')

%==========================================================================
%humidity? RH averaged over column? Or RH at source?
%==========================================================================
hres=IV_Htop_BE./fit_Htop(IV_MER_BE);
%hres=hres./fit_wshear(IV_Wshear_BE);
x=IV_RH_BE./IV_RH_ref;


subplot(2,2,3)


mask=abs(IV_latitude)<25;
plot(x(mask),hres(mask),'ks','MarkerSize',5,'MarkerFaceColor',[0.6 0 0])
hold on
mask=abs(IV_latitude)>25;
plot(x(mask),hres(mask),'k>','MarkerSize',5,'MarkerFaceColor',[0.6 0.8 1])
% plot(x,hres,'ko','MarkerSize',5,'MarkerFaceColor','k')
% plot(x(mask),hres(mask),'kp','MarkerSize',5,'MarkerFaceColor','#7a5195')
% hold on
% mask=strcmp(IV_style,'phreatomagmatic');
% plot(x(mask),hres(mask),'k^','MarkerSize',5,'MarkerFaceColor','#ffa600')
% mask=strcmp(IV_style,'phreatic' );
% plot(x(mask),hres(mask),'ko','MarkerSize',5,'MarkerFaceColor','#ef5675')
% mask=strcmp(IV_style,'Unknown' );
% plot(x(mask),hres(mask),'kd','MarkerSize',5,'MarkerFaceColor','w')
% ylim([-8 18])
% xlim([2 39])
ylabel('$\rm \overline{H}_{top}^{std}$','interpreter','latex')  
xlabel('$\rm Detrended\ relative\ humidity$','Interpreter','Latex')


mask=~isnan(hres);
[fitwind stats_RH]=fit(x(mask),hres(mask),'poly1')
[cc pp]=corrcoef(x(mask),hres(mask));cc=round(cc(1,2),2);pp(1,2)
plot(linspace(min(x),max(x),500),fitwind(linspace(min(x),max(x),500)),':','LineWidth',2,'Color',[0 0 0])
legend('Tropical events','Extra-tropical events',strcat('Linear trend: r=',num2str(cc)),'Location','NorthEast')
xlim([0.9*min(x) 1.1*max(x)])
ylim([min(hres)-0.15 max(hres)+0.15])


% for jj=1:length(eventlist3)
% ii=strcmp(IV_ivid,eventlist3{jj});text(x(ii)+2,hres(ii),eventlist3{jj},'Interpreter','none','FontSize',fontlab);
% end

% text(0.98,0.95,'Linear fit: R^2=0','Units','normalized','HorizontalAlignment','right')

%set(gca,'Yscale','log','Xscale','log')

%==========================================================================
%res plot 1
%==========================================================================
hres=IV_Htop_BE./fit_Htop(IV_MER_BE);
%hres=hres./fit_wshear(IV_Wshear_BE);
x=IV_PI_BE.*IV_W_ref./IV_N_ref;

subplot(2,2,4)
hold on
hw=plot(x(mask_weak),hres(mask_weak),'o','Color','k','MarkerFaceColor','#003f5c','MarkerSize',5)
hold on
hs=plot(x(mask_strong),hres(mask_strong),'s','Color','k','MarkerFaceColor','#bc5090','MarkerSize',5)
hu=plot(x(mask_unknown),hres(mask_unknown),'^','Color','k','MarkerFaceColor','#ffa600','MarkerSize',5)
%plotting a horizintal line highlighting a ratio of 1
max(x)
% semilogx(x,ones(size(x))*(1/1.32),'-','LineWidth',1,'Color',[0.5 0.5 0.5])
%plotting a horizintal line highlighting the mean

hold on
% ylim([-8 18])
% xlim([2 39])
ylabel('$\rm \overline{H}_{top}^{std}$','interpreter','latex')  
xlabel('$\rm Detrended\ \Pi$','Interpreter','Latex')


mask=~isnan(hres);
[fitwind stats]=fit(x(mask),hres(mask),'poly1')
[cc pp]=corrcoef(x(mask),hres(mask));cc=round(cc(1,2),2);pp(1,2)
hfit=plot(linspace(min(x),max(x),500),fitwind(linspace(min(x),max(x),500)),':','LineWidth',2,'Color',[0 0 0])
legend('Weak plumes','Strong plumes','Unknown plume morphology',strcat('Linear trend: \bf{r=',num2str(cc),'}'),'Location','NorthEast')
xlim([0.9*min(x) 1.1*max(x)])
ylim([min(hres)-0.15 max(hres)+0.15])

% 
% for jj=1:length(eventlist4)
% ii=strcmp(IV_ivid,eventlist4{jj});text(x(ii)*1.1,hres(ii),eventlist4{jj},'Interpreter','none','FontSize',fontlab);
% end

set(gca,'Xscale','log')

%==========================================================================

subplot(2,2,1)
text(-0.04,1.085,'a)','Units','normalized','FontWeight','Bold','FontSize',12)

subplot(2,2,2)
text(-0.04,1.085,'b)','Units','normalized','FontWeight','Bold','FontSize',12)

subplot(2,2,3)
text(-0.04,1.085,'c)','Units','normalized','FontWeight','Bold','FontSize',12)

subplot(2,2,4)
text(-0.04,1.085,'d)','Units','normalized','FontWeight','Bold','FontSize',12)






%==========================================================================
%res plot 1
%==========================================================================
hres=log10(IV_Htop_BE./fit_Htop(IV_MER_BE));
%hres=hres./fit_wshear(IV_Wshear_BE);
x=log10(IV_N_BE./IV_N_ref);

figure
subplot(2,2,1)

mask=abs(IV_latitude)<25;
plot(x(mask),hres(mask),'ks','MarkerSize',5,'MarkerFaceColor',[0.6 0 0])
hold on
mask=abs(IV_latitude)>25;
plot(x(mask),hres(mask),'k>','MarkerSize',5,'MarkerFaceColor',[0.6 0.8 1])

hold on
% ylim([-8 18])
% xlim([2 39])
ylabel('$\rm \overline{H}_{top}^{std}$','interpreter','latex')  
xlabel('$\rm Detrended\ Brunt$ V\"ais\"al\"a frequency','Interpreter','Latex')


mask=~isnan(hres);
[fitwind stats]=fit(x(mask),hres(mask),'poly1')
[cc pp]=corrcoef(x(mask),hres(mask));cc=round(cc(1,2),2);pp(1,2)
plot(linspace(min(x),max(x),500),fitwind(linspace(min(x),max(x),500)),':','LineWidth',2,'Color',[0 0 0])
% legend('Magmatic event','Phreatomagmatic','Phreatic','Unknown Style','Linear fit (R^2=0)','Location','NorthEast')
legend('Tropical events','Extra-tropical events',strcat('Linear trend: r=',num2str(cc)),'Location','NorthEast')
xlim([0.95*min(x) max(x)+0.05*min(x)])
ylim([min(hres)-0.15 max(hres)+0.15])



