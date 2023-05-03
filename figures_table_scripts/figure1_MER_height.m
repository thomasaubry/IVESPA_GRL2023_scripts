% 	Written by Thomas J. Aubry, May 2023.
% 	Department of Earth and Environmental Sciences, University of Exeter
%   E-mail: t.aubry@exeter.ac.uk
% 	Please cite the corresponding paper if you use this script

clear
close all

addpath('../functions/') %to access useful functions
addpath('../') %to access data
load_IVESPA %scripts reading IVESPA data


figure(1)

%==========================================================================
%Height vs MER (top height)
%==========================================================================

x=IV_TEM_BE./IV_duration_BE; % x = MER = TEM/duration
xl1=0.8*min(x);xl2=1.5*max(x); %define figure limits
dx_l=x.*((IV_TEM_UL./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5; %rough estimataion of lower and upper bound on MER using propagation rules.
dx_u=x.*((IV_TEM_UU./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;


y=IV_Htop_BE; % y = plume top height in km a.v.l.
dy=IV_Htop_U;
% filter out events that do not have the target height available.
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dx_l=dx_l(mask);dy=dy(mask);dx_u=dx_u(mask);

%fit a power law and save the model
[fit_Htop stats]=fit(x,y,'power1')
save fit_Htop fit_Htop

%define a regularly space MER grid for some predictions
xlist=logspace(log10(xl1),log10(xl2),1000);
%Calculate confidence intervals
%(https://uk.mathworks.com/help/curvefit/predint.html)
fit_int= predint(fit_Htop,xlist,0.95,'functional');
fit_int2= predint(fit_Htop,xlist,0.95,'observation');

%open figure and top-left subplot
figure(1)
subplot(3,2,1)

%create vectors with MER and height confidence interval to shade uncertainty
%using fill
X=[xlist,fliplr(xlist)];                
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];              
huncert=fill(X,Y,'r','FaceColor','#003f5c','EdgeColor','none','FaceAlpha',0.25); 
hold on

%plot the prediction uncertainty as dashed lines
huncert2=plot(xlist,fit_int2(:,1),'--','Color','#003f5c','LineWidth',1)
plot(xlist,fit_int2(:,2),'--','Color','#003f5c','LineWidth',1)

%plot the actual fit relationship, as well as the IVESPA data; put graph in
%log scale
hfit=plot(xlist,fit_Htop(xlist),'-','Color','#003f5c','LineWidth',3)
hdata=errorbar(x,y,dy,dy,dx_l,dx_u,'ko','MarkerFaceColor','k','MarkerSize',3,'CapSize',1,'LineWidth',0.03)
set(gca,'XScale','log')

%Fit the data but this time using height as the independent variable and
%using a linear fit after log transformation
[fit_Mtop stats]=fit(log10(y),log10(x),'poly1');
save fit_Mtop fit_Mtop

%Plot predicted values on a regularly spcaed height grid
ylist=linspace(0.01,100,1000);
ylist=linspace(10.^interp1(fit_Mtop(log10(ylist)),log10(ylist),log10(xl1)),10.^interp1(fit_Mtop(log10(ylist)),log10(ylist),log10(xl2)),1000);
hfitM=plot(10.^fit_Mtop(log10(ylist)),ylist,':','Color','#003f5c','LineWidth',3)

%Exctract fit coefficient values and define text showing fit equations to
%be added to legend
c_Htop = coeffvalues(fit_Htop);c_Htop=round(c_Htop,3);
c_Mtop = coeffvalues(fit_Mtop);c_Mtop=round(c_Mtop,2);
top_textH=strcat('Top height: $\rm \overline{H}_{top}$=',num2str(c_Htop(1)),'$\rm \times \overline{MER}^{',num2str(c_Htop(2)),'}$');
top_textM=strcat('Top height: log($\rm \overline{MER}$)=',num2str(c_Mtop(2)),'+',num2str(c_Mtop(1)),'$\rm \times log( \overline{H}_{top}$)');

%Set axis limits, label, subplot title and add legend
xlabel('$\rm \overline{MER} \ (kg \ s^{-1})$','Interpreter','Latex')
ylabel('$\rm \overline{H}_{top}$ (km a.v.l.)','Interpreter','Latex')
xlim([xl1 xl2])
ylim([0 42])
legend([hdata hfit huncert huncert2 hfitM],'IVESPA data','Height power law fit (H = a$\rm \times\overline{MER}^b$)','Confidence interval for height fit','Prediction interval for height fit','$\rm \overline{MER}$ log-linear fit (log($\rm \overline{MER}$) = c+d$\rm \times$log(H))','Location','Northwest','Interpreter','Latex')
title('a) Fits for top height ($\rm \overline{H}_{top}$, 130 events)','Interpreter','Latex')
%==========================================================================
%Height vs MER (spreading height)
%==========================================================================

%This section of code is exactly the same as for top height, repeated for
%spreading height
x=IV_TEM_BE./IV_duration_BE;
dx_l=x.*((IV_TEM_UL./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
dx_u=x.*((IV_TEM_UU./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
y=IV_Hspread_BE;
dy=IV_Hspread_U;
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dx_l=dx_l(mask);dy=dy(mask);dx_u=dx_u(mask);
[fit_Hspread stats]=fit(x,y,'power1')
xlist=logspace(log10(xl1),log10(xl2),1000);

fit_int= predint(fit_Hspread,xlist,0.95,'functional');
fit_int2= predint(fit_Hspread,xlist,0.95,'observation');

figure(1)
subplot(3,2,2)

X=[xlist,fliplr(xlist)];               
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];             
huncert=fill(X,Y,'r','FaceColor','#7a5195','EdgeColor','none','FaceAlpha',0.25); 
hold on
huncert2=plot(xlist,fit_int2(:,1),'--','Color','#7a5195','LineWidth',1)
plot(xlist,fit_int2(:,2),'--','Color','#7a5195','LineWidth',1)
hfit=plot(xlist,fit_Hspread(xlist),'-','Color','#7a5195','LineWidth',3)
hdata=errorbar(x,y,dy,dy,dx_l,dx_u,'ko','MarkerFaceColor','k','MarkerSize',3,'CapSize',1,'LineWidth',0.03)
set(gca,'XScale','log')


[fit_Mspread stats]=fit(log10(y),log10(x),'poly1');
ylist=linspace(0.01,100,1000);
ylist=linspace(10.^interp1(fit_Mspread(log10(ylist)),log10(ylist),log10(xl1)),10.^interp1(fit_Mspread(log10(ylist)),log10(ylist),log10(xl2)),1000);
hfitM=plot(10.^fit_Mspread(log10(ylist)),ylist,':','Color','#7a5195','LineWidth',3)
c_Hspread = coeffvalues(fit_Hspread);c_Hspread=round(c_Hspread,3);
c_Mspread = coeffvalues(fit_Mspread);c_Mspread=round(c_Mspread,2);
spread_textH=strcat('Spreading height: $\rm \overline{H}_{spr}$=',num2str(c_Hspread(1)),'$\rm \times \overline{\overline{MER}}^{',num2str(c_Hspread(2)),'}$');
spread_textM=strcat('Spreading height: log(MER)=',num2str(c_Mspread(2)),'+',num2str(c_Mspread(1)),'$\rm\times log( \overline{H}_{spr}$)');


xlabel('$\rm \overline{MER} \ (kg \ s^{-1})$','Interpreter','Latex')
ylabel('$\rm \overline{H}_{spr}$ (km a.v.l.)','Interpreter','Latex')
xlim([0.8*min(x) 1.5*max(x)])
ylim([0 42])
xlim([xl1 xl2])

title('b) Fits for spreading height ($\rm \overline{H}_{spr}$, 41 events)','Interpreter','Latex')
%==========================================================================
%Height vs MER (so2 height)
%==========================================================================
%This section of code is exactly the same as for top height, repeated for
%SO2 height

x=IV_TEM_BE./IV_duration_BE;
dx_l=x.*((IV_TEM_UL./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
dx_u=x.*((IV_TEM_UU./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
y=IV_Hso2_BE;
dy=IV_Hso2_U;
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dx_l=dx_l(mask);dy=dy(mask);dx_u=dx_u(mask);
[fit_Hso2 stats]=fit(x,y,'power1')
xlist=logspace(log10(xl1),log10(xl2),1000);

fit_int= predint(fit_Hso2,xlist,0.95,'functional');
fit_int2= predint(fit_Hso2,xlist,0.95,'observation');

figure(1)
subplot(3,2,3)

X=[xlist,fliplr(xlist)];               
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];              
huncert=fill(X,Y,'r','FaceColor','#ffa600','EdgeColor','none','FaceAlpha',0.25); 
hold on
huncert2=plot(xlist,fit_int2(:,1),'--','Color','#ffa600','LineWidth',1)
plot(xlist,fit_int2(:,2),'--','Color','#ffa600','LineWidth',1)
hfit=plot(xlist,fit_Hso2(xlist),'-','Color','#ffa600','LineWidth',3)
hdata=errorbar(x,y,dy,dy,dx_l,dx_u,'ko','MarkerFaceColor','k','MarkerSize',3,'CapSize',1,'LineWidth',0.03)
set(gca,'XScale','log')


[fit_Mso2 stats]=fit(log10(y),log10(x),'poly1');
ylist=linspace(0.01,100,1000);
ylist=linspace(10.^interp1(fit_Mso2(log10(ylist)),log10(ylist),log10(xl1)),10.^interp1(fit_Mso2(log10(ylist)),log10(ylist),log10(xl2)),1000);

hfitM=plot(10.^fit_Mso2(log10(ylist)),ylist,':','Color','#ffa600','LineWidth',3)
c_Hso2 = coeffvalues(fit_Hso2);c_Hso2=round(c_Hso2,3);
c_Mso2 = coeffvalues(fit_Mso2);c_Mso2=round(c_Mso2,2);
so2_textH=strcat('$\rm SO_2$ height: $\rm \overline{H}_{SO2}$=',num2str(c_Hso2(1)),'$\rm \times \overline{MER}^{',num2str(c_Hso2(2)),'}$');
so2_textM=strcat('$\rm SO_2$ height: log(MER)=',num2str(c_Mso2(2)),'+',num2str(c_Mso2(1)),'$\rm\times log( \overline{H}_{SO2}$)');


xlabel('$\rm \overline{MER} \ (kg \ s^{-1})$','Interpreter','Latex')
ylabel('$\rm \overline{H}_{SO2}$ (km a.v.l.)','Interpreter','Latex')
xlim([0.8*min(x) 1.5*max(x)])

ylim([0 42])
xlim([xl1 xl2])

title('c) Fits for $\rm SO_2$ height ($\rm \overline{H}_{SO2}$, 28 events)','Interpreter','Latex')
%==========================================================================
%Height vs MER (iso height)
%==========================================================================

%Load isopleth-derived height from table S2
[isodata isotext]=xlsread('TableS2.xlsx');
isoid=isotext(2:end,2); %IVESPA ID of events with isopleth height available
iso_vent=isodata(:,1)/1000; %vent height in km
iso_sampling=isodata(:,2)/1000; %sampling height in km
iso_height=isodata(:,3)+iso_sampling-iso_vent; %convert isopleth height from km above sampling level to km above vent level
iso_height_U=isodata(:,4);%uncertainty


IV_Hiso_BE=NaN(size(IV_Htop_BE));IV_Hiso_U=NaN(size(IV_Htop_BE));%create NaN-filled vectors of the same size as other IVESPA variables
for i=1:length(IV_Htop_BE)
   mask_id=strcmp(IV_ivid{i},isoid);%for each available isopleth height, find matching ID in list of IVESPA ID and fill height vectors accordingly
   if sum(mask_id)~=0
   IV_Hiso_BE(i)= iso_height(mask_id);%
   IV_Hiso_U(i)= iso_height_U(mask_id);
   end
end
%If the uncertainty is missing, take a relative uncertainty equal to the
%maximum relative uncertainty accross the dataset
IV_Hiso_U(isnan(IV_Hiso_U))=IV_Hiso_BE(isnan(IV_Hiso_U))*max(IV_Hiso_U./IV_Hiso_BE);


%The section of code below is exactly the same as for top height, repeated for
%isopleth height

x=IV_TEM_BE./IV_duration_BE;
dx_l=x.*((IV_TEM_UL./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
dx_u=x.*((IV_TEM_UU./IV_TEM_BE).^2+(IV_duration_U./IV_duration_BE).^2).^0.5;
y=IV_Hiso_BE;
dy=IV_Hiso_U;
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dx_l=dx_l(mask);dy=dy(mask);dx_u=dx_u(mask);
[fit_Hiso stats]=fit(x,y,'power1')
xlist=logspace(log10(xl1),log10(xl2),1000);
 
fit_int= predint(fit_Hiso,xlist,0.95,'functional');
fit_int2= predint(fit_Hiso,xlist,0.95,'observation');
 
figure(1)
subplot(3,2,4)

X=[xlist,fliplr(xlist)];                
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];             
huncert=fill(X,Y,'r','FaceColor','#ef5675','EdgeColor','none','FaceAlpha',0.25); 
hold on
huncert2=plot(xlist,fit_int2(:,1),'--','Color','#ef5675','LineWidth',1)
plot(xlist,fit_int2(:,2),'--','Color','#ef5675','LineWidth',1)
hfit=plot(xlist,fit_Hiso(xlist),'-','Color','#ef5675','LineWidth',3)
hdata=errorbar(x,y,dy,dy,dx_l,dx_u,'ko','MarkerFaceColor','k','MarkerSize',3,'CapSize',1,'LineWidth',0.03)
set(gca,'XScale','log')


[fit_Miso stats]=fit(log10(y),log10(x),'poly1');
ylist=linspace(0.01,100,1000);
ylist=linspace(10.^interp1(fit_Miso(log10(ylist)),log10(ylist),log10(xl1)),10.^interp1(fit_Miso(log10(ylist)),log10(ylist),log10(xl2)),1000);

hfitM=plot(10.^fit_Miso(log10(ylist)),ylist,':','Color','#ef5675','LineWidth',3)
c_Hiso = coeffvalues(fit_Hiso);c_Hiso=round(c_Hiso,3);
c_Miso = coeffvalues(fit_Miso);c_Miso=round(c_Miso,2);
iso_textH=strcat('Isopleth height: $\rm H_{iso,top}$=',num2str(c_Hiso(1)),'$\rm \times \overline{MER}^{',num2str(c_Hiso(2)),'}$');
iso_textM=strcat('Isopleth height: log($\rm \overline{MER}$)=',num2str(c_Miso(2)),'+',num2str(c_Miso(1)),'$\rm \times log(H_{iso,top}$)');


xlabel('$\rm \overline{MER} \ (kg \ s^{-1})$','Interpreter','Latex')
ylabel('$\rm H_{iso,top}$ (km a.v.l.)','Interpreter','Latex')
xlim([0.8*min(x) 1.5*max(x)])
ylim([0 42])
xlim([xl1 xl2])
title('d) Fits for isopleth height ($\rm H_{iso,top}$, 18 events)','Interpreter','Latex')

%==========================================================================
%All fit together
%==========================================================================
%In this part, I plot the power law fits obtained with MER as independent
%variables, comparing fits for all height types on the same panel


figure(1)
subplot(3,2,5)
%define a regularly spaced MER grid
xlist=logspace(log10(xl1),log10(xl2),1000);
X=[xlist,fliplr(xlist)];
%get and plot the confidence interval for the top height
fit_int= predint(fit_Htop,xlist,0.95,'functional');
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];
huncert=fill(X,Y,'r','FaceColor','#003f5c','EdgeColor','none','FaceAlpha',0.25); 
hold on
set(gca,'XScale','log')
%get and plot the confidence interval for the spreading height
fit_int= predint(fit_Hspread,xlist,0.95,'functional');
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];
huncert=fill(X,Y,'r','FaceColor','#7a5195','EdgeColor','none','FaceAlpha',0.25); 
%get and plot the confidence interval for the SO2 height
fit_int= predint(fit_Hso2,xlist,0.95,'functional');
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];
huncert=fill(X,Y,'r','FaceColor','#ffa600','EdgeColor','none','FaceAlpha',0.25); 
%get and plot the confidence interval for the isopleth height
fit_int= predint(fit_Hiso,xlist,0.95,'functional');
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];
huncert=fill(X,Y,'r','FaceColor','#ef5675','EdgeColor','none','FaceAlpha',0.25); 

%Then plot all four fit relationships (best estimate not confidence
%interval)
hiso=semilogx(xlist,fit_Hiso(xlist),'-','Color','#ef5675','LineWidth',2)
hso2=semilogx(xlist,fit_Hso2(xlist),'-','Color','#ffa600','LineWidth',2)
hspread=semilogx(xlist,fit_Hspread(xlist),'-','Color','#7a5195','LineWidth',2)
htop=semilogx(xlist,fit_Htop(xlist),'-','Color','#003f5c','LineWidth',2)

%Add the Mastin et el. 2009 + Wilson and Walker 1978 fits
hmastin=semilogx(xlist,2*(xlist/2500).^0.241,'k-.','LineWidth',1.5)
hwilson=semilogx(xlist,0.236*(xlist).^0.25,'k--','LineWidth',1.5)

%legend, axis, title, etc...
legend([htop hspread hso2 hiso hmastin hwilson],top_textH,spread_textH,so2_textH,iso_textH,'Mastin et al. (2009): H=0.304$\rm \times \overline{MER}^{0.241}$','Wilson \& Walker (1987):  H=0.236 $\rm \times \overline{MER}^{0.25}$','Location','Northwest','Interpreter','Latex')
xlabel('$\rm \overline{MER} \ (kg \ s^{-1})$','Interpreter','Latex')
ylabel('H (km a.v.l.)')
ylim([0 42])
xlim([xl1 xl2])
title('e) H-$\rm \overline{MER}$ power law fit comparison','Interpreter','Latex')



%==========================================================================
%All fit together
%==========================================================================
%This is exactly the same as the previous code section, but plotting log-
%linear fits obtained with height as independent variable
figure(1)
subplot(3,2,6)
xlist=linspace(0.1,55,1000);
X=[xlist,fliplr(xlist)];
fit_int= 10.^predint(fit_Mtop,log10(xlist),0.95,'functional');
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];
huncert=fill(X,Y,'r','FaceColor','#003f5c','EdgeColor','none','FaceAlpha',0.25); 
hold on
set(gca,'YScale','log')


fit_int= 10.^predint(fit_Mspread,log10(xlist),0.95,'functional');
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];
huncert=fill(X,Y,'r','FaceColor','#7a5195','EdgeColor','none','FaceAlpha',0.25); 

fit_int= 10.^predint(fit_Mso2,log10(xlist),0.95,'functional');
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];
huncert=fill(X,Y,'r','FaceColor','#ffa600','EdgeColor','none','FaceAlpha',0.25); 

fit_int= 10.^predint(fit_Miso,log10(xlist),0.95,'functional');
Y=[fit_int(:,1)',fliplr(fit_int(:,2)')];
huncert=fill(X,Y,'r','FaceColor','#ef5675','EdgeColor','none','FaceAlpha',0.25); 



hiso=semilogx(xlist,10.^fit_Miso(log10(xlist)),':','Color','#ef5675','LineWidth',2)
hso2=semilogx(xlist,10.^fit_Mso2(log10(xlist)),':','Color','#ffa600','LineWidth',2)
hspread=semilogx(xlist,10.^fit_Mspread(log10(xlist)),':','Color','#7a5195','LineWidth',2)
htop=semilogx(xlist,10.^fit_Mtop(log10(xlist)),':','Color','#003f5c','LineWidth',2)

hmastin=semilogx(xlist, 2500*(xlist/2).^(1/0.241),'k-.','LineWidth',1.5)
%hsparks=semilogx(xlist,1.67*(xlist/2500).^0.259,'k:','LineWidth',1)%Sparks et al. (1997)
hwilson=semilogx(xlist,(xlist/0.236).^(1/0.25),'k--','LineWidth',1.5)

ylim([xl1 2*10^9])
xlim([0 42])
ylabel('$\rm \overline{MER} \ (kg \ s^{-1})$','Interpreter','Latex')
xlabel('H (km a.v.l.)')

legend([htop hspread hso2 hiso hmastin hwilson],top_textM,spread_textM,so2_textM,iso_textM,' Mastin et al. (2009): $ \rm log(\overline{MER})$=2.15+4.15$\rm \times$log(H)','Wilson \& Walker (1987): log($\rm \overline{MER}$)=2.51+4$\rm \times$log(H)','Location','Northwest','Interpreter','Latex')
ht=title('f) $\rm \overline{MER}$-H log-linear fit comparison','Interpreter','Latex')

%==========================================================================
%Label outliers
%==========================================================================
%Label a few outliers on the top four graph. The events are selected using
%their IVESPA ID (e.g. PIN1991_02 for phase 2 of the 1991 Mt Pinatubo
%eruption)

fontlab=5;
meroffset=1.2;meroffset2=1.1
hoffset=1;
figure(1)

subplot(3,2,3)
txt='ETN2011_01';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset2,IV_Hso2_BE(ii)-hoffset,txt,'Interpreter','none','FontSize',fontlab);
txt='PIN1991_02';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)/2,IV_Hso2_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);
% txt='ELC1982_02';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset2,IV_Hso2_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);
% txt='HEL1980_02';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset2,IV_Hso2_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);
% txt='MIY2000_03';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset2,IV_Hso2_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);
txt='MER2010_01';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset2,IV_Hso2_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);
txt='NEV1985_01';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset2,IV_Hso2_BE(ii)-hoffset,txt,'Interpreter','none','FontSize',fontlab);

subplot(3,2,1)
txt='COT2015_04';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset,IV_Htop_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);
% txt='STM1902_01';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset,IV_Htop_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);
txt='BEZ1956_01';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset,IV_Htop_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);
txt='CHA2008_01';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset,IV_Htop_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);
% txt='MER2010_01';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset,IV_Htop_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);
% txt='MIY2000_03';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset,IV_Htop_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);
txt='STR2003_01';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset,IV_Htop_BE(ii)-hoffset,txt,'Interpreter','none','FontSize',fontlab);
%txt='ETN2006_01';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset,IV_Htop_BE(ii)-hoffset,txt,'Interpreter','none','FontSize',fontlab);
% txt='RED1990_06';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset,IV_Htop_BE(ii)-hoffset,txt,'Interpreter','none','FontSize',fontlab);
% 

subplot(3,2,2)
txt='POP1996_02';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii),IV_Hspread_BE(ii)+2,txt,'Interpreter','none','FontSize',fontlab);
txt='PIN1991_02';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)/2,IV_Hspread_BE(ii)+2,txt,'Interpreter','none','FontSize',fontlab);
txt='HEL1980_01';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)/2,IV_Hspread_BE(ii)+2,txt,'Interpreter','none','FontSize',fontlab);
% txt='ELC1982_03';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)/4,IV_Hspread_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);
% txt='ELC1982_02';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset,IV_Hspread_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);


subplot(3,2,4)
txt='PIN1991_02';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)/2,IV_Hiso_BE(ii)-hoffset,txt,'Interpreter','none','FontSize',fontlab);
txt='ETN2011_01';ii=strcmp(IV_ivid,txt);text(IV_MER_BE(ii)*meroffset,IV_Hiso_BE(ii)+hoffset,txt,'Interpreter','none','FontSize',fontlab);


subplot(3,2,6)
set(gca,'YTick',10.^[2:2:8])


