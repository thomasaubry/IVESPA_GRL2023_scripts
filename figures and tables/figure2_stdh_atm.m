% 	Written by Thomas J. Aubry, May 2023.
% 	Department of Earth and Environmental Sciences, University of Exeter
%   E-mail: t.aubry@exeter.ac.uk
% 	Please cite the corresponding paper if you use this script

clear
close all

addpath('../functions/') %to access useful functions
addpath('../') %to access data
load_IVESPA %scripts reading IVESPA data
load('fit_Htop.mat') %load top height fit from figure 1



%define lists of outliers to label on figure 2 panels
eventlist={'RED1990_06','STR2003_01','RED1990_07','REV2002_01','NEV1985_01','RED1990_02','QUI1932_01','STM1902_01','ETN2013_01','ANA2003_02','EYJ2010_04','CHA2008_03','SHM1997_01','MER2010_01','CHA2008_01','BEZ1956_01'};
eventlist1={'MER2010_01','CHA2008_01','BEZ1956_01','MIY2000_01','COT2015_01'};
eventlist2={'MER2010_01','CHA2008_01','MIY2000_01','COT2015_01','RED1990_05'};
eventlist3={'MER2010_01','CHA2008_01','MIY2000_01','COT2015_01','GRI2004_02','NEV1985_01'};
eventlist4={'MER2010_01','CHA2008_01','MIY2000_01','COT2015_01','ETN2016_02'};

%define a few plotting parameters
fontlab=5;
meroffset=1.2;meroffset2=1.1
hoffset=1;

%==========================================================================
%Standardized height vs wind speed
%==========================================================================
hres=IV_Htop_BE./fit_Htop(IV_MER_BE);%define hres as the standardized top height
x=IV_W_BE; % define x as wind speed

figure
subplot(2,2,2)
mask=abs(IV_latitude)<25;%first plot tropical events
plot(x(mask),hres(mask),'ks','MarkerSize',5,'MarkerFaceColor',[0.6 0 0])
hold on
mask=abs(IV_latitude)>25;%then high-latitude events
plot(x(mask),hres(mask),'k>','MarkerSize',5,'MarkerFaceColor',[0.6 0.8 1])
hold on




mask=~isnan(hres);%create a mask to remove event with no top height estimate
[fitwind stats]=fit(x(mask),hres(mask),'poly1');%linear fit of standardized height vs wind
[cc pp]=corrcoef(x(mask),hres(mask));cc=round(cc(1,2),2);pp(1,2)% correlation coefficient for standardized height vs wind
%add linear trend to graph
plot(linspace(min(x),max(x),500),fitwind(linspace(min(x),max(x),500)),':','LineWidth',2,'Color',[0 0 0])

%legend and axis limit and labels
legend('Tropical events','Extra-tropical events',strcat('Linear trend: r=',num2str(cc)),'Location','NorthEast')
xlim([0.9*min(x) 1.1*max(x)])
ylim([min(hres)-0.15 max(hres)+0.15])
ylabel('$\rm \overline{H}_{top}^{std}$','interpreter','latex')   
xlabel('Wind speed $\rm \overline{W}$ ($\rm m \ s^{-1}$)','Interpreter','Latex')


%label a few outliers
for jj=1:length(eventlist2)
ii=strcmp(IV_ivid,eventlist2{jj});text(x(ii)+1,hres(ii),eventlist2{jj},'Interpreter','none','FontSize',fontlab);
end


%==========================================================================
%Standardized height vs Brunt Vaisala
%==========================================================================
%code section not commented but same as for wind speed

hres=IV_Htop_BE./fit_Htop(IV_MER_BE);
x=IV_N_BE;

subplot(2,2,1)

mask=abs(IV_latitude)<25;
plot(x(mask),hres(mask),'ks','MarkerSize',5,'MarkerFaceColor',[0.6 0 0])
hold on
mask=abs(IV_latitude)>25;
plot(x(mask),hres(mask),'k>','MarkerSize',5,'MarkerFaceColor',[0.6 0.8 1])

hold on
ylabel('$\rm \overline{H}_{top}^{std}$','interpreter','latex')   
xlabel('Brunt V\"ais\"al\"a frequency $\rm \overline{N}$ ($\rm s^{-1}$)','Interpreter','Latex')


mask=~isnan(hres);
[fitwind stats]=fit(x(mask),hres(mask),'poly1')
[cc pp]=corrcoef(x(mask),hres(mask));cc=round(cc(1,2),2);pp(1,2)
plot(linspace(min(x),max(x),500),fitwind(linspace(min(x),max(x),500)),':','LineWidth',2,'Color',[0 0 0])
legend('Tropical events','Extra-tropical events',strcat('Linear trend: r=',num2str(cc)),'Location','NorthEast')
xlim([0.009 0.021])
ylim([min(hres)-0.15 max(hres)+0.15])



for jj=1:length(eventlist1)
ii=strcmp(IV_ivid,eventlist1{jj});text(x(ii)+0.0003,hres(ii),eventlist1{jj},'Interpreter','none','FontSize',fontlab);
end
%set(gca,'Yscale','log','Xscale','log')

%==========================================================================
%Standardized height vs relative humidity
%==========================================================================
%code section not commented but same as for wind speed

hres=IV_Htop_BE./fit_Htop(IV_MER_BE);
x=IV_RH_BE;

subplot(2,2,3)

mask=abs(IV_latitude)<25;
plot(x(mask),hres(mask),'ks','MarkerSize',5,'MarkerFaceColor',[0.6 0 0])
hold on
mask=abs(IV_latitude)>25;
plot(x(mask),hres(mask),'k>','MarkerSize',5,'MarkerFaceColor',[0.6 0.8 1])
ylabel('$\rm \overline{H}_{top}^{std}$','interpreter','latex')   
xlabel('Relative humidity ($\rm \%$)','Interpreter','Latex')


mask=~isnan(hres);
[fitwind stats_RH]=fit(x(mask),hres(mask),'poly1')
[cc pp]=corrcoef(x(mask),hres(mask));cc=round(cc(1,2),2);pp(1,2)
plot(linspace(min(x),max(x),500),fitwind(linspace(min(x),max(x),500)),':','LineWidth',2,'Color',[0 0 0])
legend('Tropical events','Extra-tropical events',strcat('Linear trend: r=',num2str(cc)),'Location','NorthEast')
xlim([0.9*min(x) 1.1*max(x)])
ylim([min(hres)-0.15 max(hres)+0.15])


for jj=1:length(eventlist3)
ii=strcmp(IV_ivid,eventlist3{jj});text(x(ii)+2,hres(ii),eventlist3{jj},'Interpreter','none','FontSize',fontlab);
end


%==========================================================================
%Standardized height vs PI
%==========================================================================

z=IV_morpho(mask); %load plume morphology (weak/strong)
mask_weak=strcmp('weak',z);%create a mask for weak plumes
mask_strong=strcmp('strong',z);%create a mask for strong plumes
mask_unknown=~mask_weak & ~mask_strong;%create a mask for unknown plumes


%code section below same as for wind speed, except that events are
%color-coded according to plume morphology instead of eruption latitude

hres=IV_Htop_BE./fit_Htop(IV_MER_BE);
x=IV_PI_BE;

subplot(2,2,4)
hold on
hw=plot(x(mask_weak),hres(mask_weak),'o','Color','k','MarkerFaceColor','#003f5c','MarkerSize',5)
hold on
hs=plot(x(mask_strong),hres(mask_strong),'s','Color','k','MarkerFaceColor','#bc5090','MarkerSize',5)
hu=plot(x(mask_unknown),hres(mask_unknown),'^','Color','k','MarkerFaceColor','#ffa600','MarkerSize',5)

max(x)


hold on

ylabel('$\rm \overline{H}_{top}^{std}$','interpreter','latex')   
xlabel('Wind entrainment/plume rise timescale ratio $\rm \Pi$','Interpreter','Latex')


mask=~isnan(hres);
[fitwind stats]=fit(x(mask),hres(mask),'poly1')
[cc pp]=corrcoef(x(mask),hres(mask));cc=round(cc(1,2),2);pp(1,2)
hfit=plot(linspace(min(x),max(x),500),fitwind(linspace(min(x),max(x),500)),':','LineWidth',2,'Color',[0 0 0])
legend('Weak plumes','Strong plumes','Unknown plume morphology',strcat('Linear trend: \bf{r=',num2str(cc),'}'),'Location','NorthEast')
xlim([0.9*min(x) 1.1*max(x)])
ylim([min(hres)-0.15 max(hres)+0.15])


for jj=1:length(eventlist4)
ii=strcmp(IV_ivid,eventlist4{jj});text(x(ii)*1.1,hres(ii),eventlist4{jj},'Interpreter','none','FontSize',fontlab);
end

set(gca,'Xscale','log')

%==========================================================================
%add panel number
subplot(2,2,1)
text(-0.04,1.085,'a)','Units','normalized','FontWeight','Bold','FontSize',12)

subplot(2,2,2)
text(-0.04,1.085,'b)','Units','normalized','FontWeight','Bold','FontSize',12)

subplot(2,2,3)
text(-0.04,1.085,'c)','Units','normalized','FontWeight','Bold','FontSize',12)

subplot(2,2,4)
text(-0.04,1.085,'d)','Units','normalized','FontWeight','Bold','FontSize',12)



