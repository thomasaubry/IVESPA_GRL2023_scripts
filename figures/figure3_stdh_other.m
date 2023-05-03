% 	Written by Thomas J. Aubry, May 2023.
% 	Department of Earth and Environmental Sciences, University of Exeter
%   E-mail: t.aubry@exeter.ac.uk
% 	Please cite the corresponding paper if you use this script

clear
close all

addpath('../functions/') %to access useful functions
addpath('../') %to access data
load_IVESPA %scripts reading IVESPA data
load('fit_Htop.mat');%load top height fit from figure 1
fit_Htop_ref=fit_Htop;

%==========================================================================
%Define a few variables
%==========================================================================
%MER and its uncertainty
x=IV_MER_BE;
dx_l=x.*((IV_TEM_UL./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
dx_u=x.*((IV_TEM_UU./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;

%top height and its uncertainty
y=IV_Htop_BE;
dy=IV_Htop_U;

%mask all variables to remove events with no top height estimates
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dx_l=dx_l(mask);dy=dy(mask);dx_u=dx_u(mask);
meanwind=IV_W_BE(mask);
IV_heighttech=IV_heighttech(mask);
IV_ivid=IV_ivid(mask);

%create a regularly-spaced list of MER
xlist=logspace(min(log10(x)),max(log10(x)),1000);


fontlab=5;%define font size


%==========================================================================
%Standardized height distribution for each eruption region
%==========================================================================
hres=y./fit_Htop_ref(x);%define hres as standardized top height
reg=IV_region(mask);%load and mask eruption regions

%define a unique region list and sort it by median standardized height
%value
reglist=unique(reg);
reglist=reglist([2 4 8 1 7 10 5 3 6 9]);

figure
subplot(4,2,[1 2])

for i=1:length(reglist)
    
    %for each region, find all corresponding events and select their wind
    %speeds and standardized heights
    mask=strcmp(reg,reglist{i});
    xv=meanwind(mask);
    yv=hres(mask);
    
    %make a box plot of the standardized height distribution for this
    %region
    boxchart(i*ones(size(yv)),yv,'BoxFaceColor','#003f5c','MarkerStyle','none')
    
    %also plot actual datapoints to visualize the distribution. Add small
    %random number to x value to vizualize more easily distinct data points
    hold on
    plot(i+(rand(size(yv))-0.5)/10,yv,'.','MarkerSize',15,'Color','#003f5c')
    
    %label the number of events in the region (n)
    text(i+0.1,2.9,strcat('n=',num2str(length(yv))),'interpreter','latex')
    
    %label the p-value (p) from a Mann-Whitney U test testing if the median
    %standarzed height from the region differs from that of all other
    %region
    pmwu=round(ranksum(hres(mask),hres(~mask)),2);
    if pmwu<=0.05 %label in bold if significant at 95% level
        text(i+0.1,2.6,strcat('\textbf{$\rm p$=',num2str(pmwu),'}'),'fontweight','bold','interpreter','latex')
    else
        text(i+0.1,2.6,strcat('$\rm p$=',num2str(pmwu)),'interpreter','latex')
    end
    
    %calculate correlation coefficient between standardized height and wind
    %speed accross region; calculated on a log scale as dependence expected
    %to be non-linear
    
    [aa pp]=corrcoef(log10(xv),log10(yv));
    if pp(2,1)<=0.05 %label correlation coefficient in bold if significant at 95% level
        text(i+0.1,2.3,strcat('\textbf{$\rm r$=',num2str(round(aa(2,1),2)),'}'),'fontweight','bold','interpreter','latex')
    else
        text(i+0.1,2.3,strcat('$\rm r$=',num2str(round(aa(2,1),2))),'interpreter','latex')
    end
    
    nn(i)=length(yv);
end

%plot the standardized height = 1 line
plot([0 12],[1 1],'k--')

%shorten a few region labels
reglist{1}='Ctrl. Amer.';
reglist{5}='Kamch. & Aleut.';
reglist{6}='SE Asia';

%axis limit and labels, and title
ylim([0.2 3.2])
xlim([0.5 length(reglist)+0.5])
set(gca,'XTick',1:length(reglist),'XTickLabel',reglist)
xtickangle(15)
ylabel('$\rm \overline{H}_{top}^{std}$','interpreter','latex')
title('\bf a) Influence of volcanic region','interpreter','latex','FontSize',11)





%==========================================================================
%Standardized height distribution for each height measurement technique
%==========================================================================
subplot(4,2,[3 4])
hres=y./fit_Htop_ref(x);%define hres as standardized top height

%define height measurement technique as tech, then do a bit of cleaning as
%combinations of techniques are commonly permutated when entered in IVESPA,
%e.g. satelitte and visual might be noted 's,v' or 'v,s'
tech=IV_heighttech;
tech(strcmp(tech,'v,g,s'))={'s,g,v'};tech(strcmp(tech,'g,v,s'))={'s,g,v'};tech(strcmp(tech,'g,s,v'))={'s,g,v'};
tech(strcmp(tech,'s,v'))={'v,s'};
tech(strcmp(tech,'v,g'))={'g,v'};
tech(strcmp(tech,'u,s'))={'s,u'};
tech(strcmp(tech,'r'))={'g'};
tech(strcmp(tech,'r,s'))={'g,s'};
tech(strcmp(tech,'v,u'))={'v'};
tech(strcmp(tech,'s,u'))={'s'};
tech(strcmp(tech,'g,u'))={'g'};
tech(strcmp(tech,'g,s'))={'s,g'};

%define a unique list of technique and sort it by median standardized
%height value
techlist=unique(tech);
techlist=techlist([4 3 5 2 1 6 7 8]);


%the code below is the exact same as for the region plot
for i=1:length(techlist)
    
    mask=strcmp(tech,techlist{i});
    xv=meanwind(mask);
    yv=hres(mask);
    boxchart(i*ones(size(yv)),yv,'BoxFaceColor','#003f5c','MarkerStyle','none')
    
    hold on
    plot(i+(rand(size(yv))-0.5)/10,yv,'.','MarkerSize',15,'Color','#003f5c')
    text(i+0.1,2.9,strcat('n=',num2str(length(yv))),'interpreter','latex')
    
    pmwu=round(ranksum(hres(mask),hres(~mask)),2);
    if pmwu<=0.05
        text(i+0.1,2.6,strcat('\textbf{$\rm p$=',num2str(pmwu),'}'),'fontweight','bold','interpreter','latex')
    else
        text(i+0.1,2.6,strcat('$\rm p$=',num2str(pmwu)),'interpreter','latex')
    end
    
    
    [aa pp]=corrcoef(log10(xv),log10(yv));

    if pp(2,1)<=0.1
        text(i+0.1,2.3,strcat('\textbf{$\rm r$=',num2str(round(aa(2,1),2)),'}'),'fontweight','bold','interpreter','latex')
    else
        text(i+0.1,2.3,strcat('$\rm r$=',num2str(round(aa(2,1),2))),'interpreter','latex')
    end
    
    
    nn(i)=length(yv);
end
plot([0.5 length(techlist)+0.5],[1 1],'k--')


techlist{2}='s (satellite)';
techlist{5}='g (ground)';
techlist{6}='unknown';
techlist{7}='v (visual)';
ylim([0.2 3.2])
xlim([0.5 length(techlist)+0.5])
set(gca,'XTick',1:length(techlist),'XTickLabel',techlist)
xtickangle(15)
ylabel('$\rm \overline{H}_{top}^{std}$','interpreter','latex')
title('\bf b) Influence of $\rm \overline{H}_{top}$ measurement technique','interpreter','latex','FontSize',11)


%==========================================================================
%MER-top height power law fir for specific height measurement techniques
%==========================================================================

subplot(4,2,[6 8])

%find all events for which top height was measured using satellite-only or
%a combination of satelitte and ground-based techniques
mask=strcmp(tech,'s,g') | strcmp(tech,'s') | strcmp(tech,'s,g,v');
%find best MER-top height power law fit for these events and calculate
%confidence interval
[fit_Htop_s stats]=fit(x(mask),y(mask),'power1')
fit_int= predint(fit_Htop_s,xlist,0.95,'functional');

%plot confidence interval first
X=[xlist,fliplr(xlist)];                
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];
huncert=fill(X,Y,'r','FaceColor','#003f5c','EdgeColor','none','FaceAlpha',0.25);
hold on

%repeat the above for top height measured using visual measurements
mask=strcmp(tech,'v') | strcmp(tech,'v,s');
[fit_Htop stats]=fit(x(mask),y(mask),'power1')
fit_int= predint(fit_Htop,xlist,0.95,'functional');
X=[xlist,fliplr(xlist)];                
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];             
huncert=fill(X,Y,'r','FaceColor','#ffa600','EdgeColor','none','FaceAlpha',0.25);

%now plot actual fit line, including the fit previously obtained with all
%IVESPA events with a top height.
hfit_a=plot(xlist,fit_Htop_ref(xlist),':','Color','k','LineWidth',2)
hold on
hfit_s=plot(xlist,fit_Htop_s(xlist),'-','Color','#003f5c','LineWidth',2)
hold on
hfit_v=plot(xlist,fit_Htop(xlist),'-','Color','#ffa600','LineWidth',2)
hold on

%Next plot the actual datapoint for satelite measurement
mask=strcmp(tech,'s,g') | strcmp(tech,'s') | strcmp(tech,'s,g,v');
hdata_s=errorbar(x(mask),y(mask),dy(mask),dy(mask),dx_l(mask),dx_u(mask),'ko','MarkerFaceColor','#003f5c','MarkerSize',3,'CapSize',1,'LineWidth',0.03)
%also retrieve fit coefficients and define text to be written in legend
c_Htop_s = coeffvalues(fit_Htop_s);c_Htop_s=round(c_Htop_s,3);
txt_sat=strcat('$\rm \overline{H}_{top}=',num2str(c_Htop_s(1)),'\times \overline{MER}^{',num2str(c_Htop_s(2)),'}$');

%repeat the above for the visual measurements
mask=strcmp(tech,'v') | strcmp(tech,'v,s');
hdata_v=errorbar(x(mask),y(mask),dy(mask),dy(mask),dx_l(mask),dx_u(mask),'k^','MarkerFaceColor','#ffa600','MarkerSize',3,'CapSize',1,'LineWidth',0.03)
c_Htop = coeffvalues(fit_Htop);c_Htop=round(c_Htop,3);
txt_vis=strcat('$\rm \overline{H}_{top}=',num2str(c_Htop(1)),'\times \overline{MER}^{',num2str(c_Htop(2)),'}$');

%add legend to plot
legend([hfit_a hfit_s hdata_s hfit_v hdata_v],'Fit for all events','Satellite or satellite+ground',strcat('events and fit:',txt_sat),'Visual or visual+satellite',strcat('events and fit:',txt_vis),'interpreter','latex')

%label a few outliers using their IVESPA ID
eventlist={'MER2010_01','CHA2008_01','BEZ1956_01','COT2015_01','STR2003_01','CNG1999_01','MIY2000_01'};
for jj=1:length(eventlist)
    ii=strcmp(IV_ivid,eventlist{jj});text(x(ii)*1.1,y(ii)+1,eventlist{jj},'Interpreter','none','FontSize',fontlab);
end

%set axis limits/labels and title
set(gca,'XScale','log')
xlabel('$\rm \overline{MER}$ (kg $\rm s^{-1}$)','interpreter','latex')
ylabel('$\rm \overline{H}_{top}$ (km a.v.l.)','interpreter','latex')
xlim([0.5*min(x) 2*max(x)])
ylim([0 max(y)*1.3])
title('\bf d) Fits for specific height measurement techniques','interpreter','latex','FontSize',11)


%==========================================================================
%MER-top height power law fir for specific regions
%==========================================================================
%Exact same as for previous plot, but repeating for specific regions
%instead of specific height measurement techniques
subplot(4,2,[5 7])
mask=strcmp(reg,'Central America');
[fit_Htop_ca stats]=fit(x(mask),y(mask),'power1')

fit_int= predint(fit_Htop_ca,xlist,0.95,'functional');
X=[xlist,fliplr(xlist)];                
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];              
huncert=fill(X,Y,'r','FaceColor','#003f5c','EdgeColor','none','FaceAlpha',0.25);
hold on

mask=strcmp(reg,'Kamchatcka & Aleutians');
[fit_Htop_kam stats]=fit(x(mask),y(mask),'power1')
fit_int= predint(fit_Htop_kam,xlist,0.95,'functional');
X=[xlist,fliplr(xlist)];                
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];              
huncert=fill(X,Y,'r','FaceColor','#bc5090','EdgeColor','none','FaceAlpha',0.25);
hold on


mask=strcmp(reg,'Redoubt');
[fit_Htop_red stats]=fit(x(mask),y(mask),'power1')
fit_int= predint(fit_Htop_red,xlist,0.95,'functional');
X=[xlist,fliplr(xlist)];                
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];              
huncert=fill(X,Y,'r','FaceColor','#ffa600','EdgeColor','none','FaceAlpha',0.25);
hold on

mask=strcmp(reg,'Central America');

hfit_all=plot(xlist,fit_Htop_ref(xlist),':','Color','k','LineWidth',2)
hold on
hfit_ca=plot(xlist,fit_Htop_ca(xlist),'-','Color','#003f5c','LineWidth',2)
hold on
hdata_ca=errorbar(x(mask),y(mask),dy(mask),dy(mask),dx_l(mask),dx_u(mask),'ko','MarkerFaceColor','#003f5c','MarkerSize',3,'CapSize',1,'LineWidth',0.03)
c_Htop_ca = coeffvalues(fit_Htop_ca);c_Htop_ca=round(c_Htop_ca,3);
txt_amer=strcat('$\rm \overline{H}_{top}=',num2str(c_Htop_ca(1)),'\times \overline{MER}^{',num2str(c_Htop_ca(2)),'}$');




mask=strcmp(reg,'Kamchatcka & Aleutians');
[fit_Htop_kam stats]=fit(x(mask),y(mask),'power1')

hfit_kam=plot(xlist,fit_Htop_kam(xlist),'-','Color','#bc5090','LineWidth',2)
hold on
hdata_kam=errorbar(x(mask),y(mask),dy(mask),dy(mask),dx_l(mask),dx_u(mask),'ks','MarkerFaceColor','#bc5090','MarkerSize',3,'CapSize',1,'LineWidth',0.03)
c_Htop_kam = coeffvalues(fit_Htop_kam);c_Htop_kam=round(c_Htop_kam,3);
txt_kam=strcat('$\rm \overline{H}_{top}=',num2str(c_Htop_kam(1)),'\times \overline{MER}^{',num2str(c_Htop_kam(2)),'}$');



mask=strcmp(reg,'Redoubt');
[fit_Htop_red stats]=fit(x(mask),y(mask),'power1')
hfit_red=plot(xlist,fit_Htop_red(xlist),'-','Color','#ffa600','LineWidth',2)
hold on
c_Htop_red = coeffvalues(fit_Htop_red);c_Htop_red=round(c_Htop_red,3);
txt_redoubt=strcat('$\rm \overline{H}_{top}=',num2str(c_Htop_red(1)),'\times \overline{MER}^{',num2str(c_Htop_red(2)),'}$');

hdata_red=errorbar(x(mask),y(mask),dy(mask),dy(mask),dx_l(mask),dx_u(mask),'k^','MarkerFaceColor','#ffa600','MarkerSize',3,'CapSize',1,'LineWidth',0.03)

legend([hfit_all hfit_ca hdata_ca hfit_kam hdata_kam hfit_red hdata_red],'Fit for all events','Central America events',strcat('and fit:',txt_amer),'Kamchatka \& Aleutians events',strcat('and fit:',txt_kam),'Redoubt events and fit:',txt_redoubt,'interpreter','latex')

set(gca,'XScale','log')

xlabel('$\rm \overline{MER}$ (kg $\rm s^{-1}$)','interpreter','latex')
ylabel('$\rm \overline{H}_{top}$ (km a.v.l.)','interpreter','latex')
xlim([0.5*min(x) 2*max(x)])
ylim([0 max(y)*1.3])

title('\bf c) Fits for specific regions','interpreter','latex','FontSize',11)

eventlist={'RED1990_06','NEV1985_01','SHM1997_01','BEZ1956_01','COT2015_01'};
for jj=1:length(eventlist)
    ii=strcmp(IV_ivid,eventlist{jj});text(x(ii)*1.1,y(ii)+1,eventlist{jj},'Interpreter','none','FontSize',fontlab);
end


