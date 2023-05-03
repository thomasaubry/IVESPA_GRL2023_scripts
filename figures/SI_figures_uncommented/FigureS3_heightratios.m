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
%Top-spreading ratio vs MER and morphology
%==========================================================================

%For convenience I define variables x=MER and y=top height, so that I can
%later copy paste this whole bloc of code and just change the definition of
%x and y to other variables

%x=MER=TEM/duration
x=IV_PI_BE;
dx=IV_PI_U;
%y=ratio of spreading height to top height (all avl)
y=(IV_Hspread_BE)./(IV_Htop_BE);
%here I propagate the uncertainty on y
dy=y.*((IV_Htop_U./IV_Htop_BE).^2+(IV_Hspread_U./IV_Hspread_BE).^2).^0.5;
%below I remove all NaN values (i.e. event for which we don't have both the
%top and spreading height)
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dy=dy(mask);dx=dx(mask);
z=IV_morpho(mask);
sum(~isnan(y))
%and then I create some mask for weak/strong/unknown plume to plot them
%with different color code
mask_weak=strcmp('weak',z);
mask_strong=strcmp('strong',z);
mask_unknown=~mask_weak & ~mask_strong;


%everything below is plotting stuff...
figure(3)
%I devide the figure in 3 panel and will use the first one for the
%spreading to top height ratio.
subplot(1,3,1)

errorbar(x,y,dy,dy,dx,dx,'k.','MarkerSize',0.1,'CapSize',2,'LineWidth',0.2)
set(gca,'XScale','log')
hold on
hw=semilogx(x(mask_weak),y(mask_weak),'o','Color','k','MarkerFaceColor','#003f5c','MarkerSize',6)
hold on
hs=semilogx(x(mask_strong),y(mask_strong),'s','Color','k','MarkerFaceColor','#bc5090','MarkerSize',6)
hu=semilogx(x(mask_unknown),y(mask_unknown),'^','Color','k','MarkerFaceColor','#ffa600','MarkerSize',6)
%plotting a horizintal line highlighting a ratio of 1
hsl=semilogx([0.5*min(x) 2*max(x)],[1 1],':','LineWidth',2,'Color',[0.25 0.25 0.25])
% semilogx(x,ones(size(x))*(1/1.32),'-','LineWidth',1,'Color',[0.5 0.5 0.5])
%plotting a horizintal line highlighting the mean
hml=semilogx([0.5*min(x) 2*max(x)],[mean(y) mean(y)],'k-','LineWidth',3)

%just adding legend/axis labels/etc
legend([hw,hs,hu,hml,hsl],'Weak plume','Strong plume','Unknown plume morphology','mean ratio value','y=1 line')

xlabel({'$\rm \Pi$ (wind entrainment/','plume rise timescale ratio)'},'Interpreter','Latex')
ylabel('$\rm \overline{H}_{spr}/\overline{H}_{top}$','Interpreter','Latex')
xlim([0.5*min(x) 2*max(x)])
ylim([0 2.5])
meanval(1)=mean(y);

% title('a) $\rm \overline{H}_{spr}/\overline{H}_{top}$ ratio','Interpreter','Latex')
%After that, I essentially copied past the same code to do the same for the
%SO2/top ratio and isopleth/top ratio
%==========================================================================
%Top-SO2 ratio vs MER and morphology
%==========================================================================

%x=MER=TEM/duration
x=IV_PI_BE;
dx=IV_PI_U;
%y=ratio of spreading height to top height (all avl)
y=(IV_Hso2_BE)./(IV_Htop_BE);
%here I propagate the uncertainty on y
dy=y.*((IV_Hso2_U./IV_Htop_BE).^2+(IV_Hso2_U./IV_Hso2_BE).^2).^0.5;
%below I remove all NaN values (i.e. event for which we don't have both the
%top and spreading height)
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dy=dy(mask);dx=dx(mask);
z=IV_morpho(mask);
sum(~isnan(y))
%and then I create some mask for weak/strong/unknown plume to plot them
%with different color code
mask_weak=strcmp('weak',z);
mask_strong=strcmp('strong',z);
mask_unknown=~mask_weak & ~mask_strong;


%everything below is plotting stuff...
figure(3)
%I devide the figure in 3 panel and will use the first one for the
%spreading to top height ratio.
subplot(1,3,2)

errorbar(x,y,dy,dy,dx,dx,'k.','MarkerSize',0.1,'CapSize',2,'LineWidth',0.2)
set(gca,'XScale','log')
hold on
hw=semilogx(x(mask_weak),y(mask_weak),'o','Color','k','MarkerFaceColor','#003f5c','MarkerSize',6)
hold on
hs=semilogx(x(mask_strong),y(mask_strong),'s','Color','k','MarkerFaceColor','#bc5090','MarkerSize',6)
hu=semilogx(x(mask_unknown),y(mask_unknown),'^','Color','k','MarkerFaceColor','#ffa600','MarkerSize',6)
%plotting a horizintal line highlighting a ratio of 1
hsl=semilogx([0.5*min(x) 2*max(x)],[1 1],':','LineWidth',2,'Color',[0.25 0.25 0.25])
% semilogx(x,ones(size(x))*(1/1.32),'-','LineWidth',1,'Color',[0.5 0.5 0.5])
%plotting a horizintal line highlighting the mean
hml=semilogx([0.5*min(x) 2*max(x)],[mean(y) mean(y)],'k-','LineWidth',3)

%just adding legend/axis labels/etc
%legend([hw,hs,hu,hml,hsl],'Weak plume','Strong plume','Unknown plume morphology','mean ratio value','y=1 line')

xlabel({'$\rm \Pi$ (wind entrainment/','plume rise timescale ratio)'},'Interpreter','Latex')
ylabel('$\rm \overline{H}_{SO2}/\overline{H}_{top}$','Interpreter','Latex')
xlim([0.5*min(x) 2*max(x)])
ylim([0 2.5])

% title('b) $\rm \overline{H}_{SO2}/\overline{H}_{top}$ ratio','Interpreter','Latex')

meanval(2)=mean(y);
%==========================================================================
%Top-isopleth ratio vs MER and morphology
%==========================================================================
[isodata isotext]=xlsread('TableS2.xlsx');
isoid=isotext(2:end,2);
iso_vent=isodata(:,1)/1000;
iso_sampling=isodata(:,2)/1000;
iso_height=isodata(:,3)+iso_sampling-iso_vent;
iso_height_U=isodata(:,4);

IV_Hiso_BE=NaN(size(IV_Htop_BE));IV_Hiso_U=NaN(size(IV_Htop_BE));
for i=1:length(IV_Htop_BE)
   mask_id=strcmp(IV_ivid{i},isoid);
   if sum(mask_id)~=0
   IV_Hiso_BE(i)= iso_height(mask_id);
   IV_Hiso_U(i)= iso_height_U(mask_id);
   end
end
IV_Hiso_U(isnan(IV_Hiso_U))=IV_Hiso_BE(isnan(IV_Hiso_U))*max(IV_Hiso_U./IV_Hiso_BE);

%x=MER=TEM/duration
x=IV_PI_BE;
dx=IV_PI_U;
%y=ratio of spreading height to top height (all avl)
y=(IV_Hiso_BE)./(IV_Htop_BE);
sum(~isnan(y))
%here I propagate the uncertainty on y
dy=y.*((IV_Htop_U./IV_Htop_BE).^2+(IV_Hiso_U./IV_Hiso_BE).^2).^0.5;
%below I remove all NaN values (i.e. event for which we don't have both the
%top and spreading height)
mask=~isnan(x) & ~isnan(y);
x=x(mask);y=y(mask);dy=dy(mask);dx=dx(mask);
z=IV_morpho(mask);

%and then I create some mask for weak/strong/unknown plume to plot them
%with different color code
mask_weak=strcmp('weak',z);
mask_strong=strcmp('strong',z);
mask_unknown=~mask_weak & ~mask_strong;


%everything below is plotting stuff...
figure(3)
%I devide the figure in 3 panel and will use the first one for the
%spreading to top height ratio.
subplot(1,3,3)

errorbar(x,y,dy,dy,dx,dx,'k.','MarkerSize',0.1,'CapSize',2,'LineWidth',0.2)
set(gca,'XScale','log')
hold on
hw=semilogx(x(mask_weak),y(mask_weak),'o','Color','k','MarkerFaceColor','#003f5c','MarkerSize',6)
hold on
hs=semilogx(x(mask_strong),y(mask_strong),'s','Color','k','MarkerFaceColor','#bc5090','MarkerSize',6)
hu=semilogx(x(mask_unknown),y(mask_unknown),'^','Color','k','MarkerFaceColor','#ffa600','MarkerSize',6)
%plotting a horizintal line highlighting a ratio of 1
hsl=semilogx([0.5*min(x) 2*max(x)],[1 1],':','LineWidth',2,'Color',[0.25 0.25 0.25])
% semilogx(x,ones(size(x))*(1/1.32),'-','LineWidth',1,'Color',[0.5 0.5 0.5])
%plotting a horizintal line highlighting the mean
hml=semilogx([0.5*min(x) 2*max(x)],[mean(y) mean(y)],'k-','LineWidth',3)

%just adding legend/axis labels/etc
%legend([hw,hs,hu,hml,hsl],'Weak plume','Strong plume','Unknown plume morphology','mean ratio value','y=1 line')

xlabel({'$\rm \Pi$ (wind entrainment/','plume rise timescale ratio)'},'Interpreter','Latex')
ylabel('$\rm H_{iso,top}/\overline{H}_{top}$','Interpreter','Latex')
xlim([0.5*min(x) 2*max(x)])
ylim([0 2.5])

meanval(3)=mean(y);

% title('c) $\rm H_{iso,top}/\overline{H}_{top}$ ratio','Interpreter','Latex')

subplot(1,3,1)
text(-0.2,1.01,'a)','Units','normalized','FontWeight','Bold','FontSize',11)

subplot(1,3,2)
text(-0.2,1.01,'b)','Units','normalized','FontWeight','Bold','FontSize',11)

subplot(1,3,3)
text(-0.2,1.01,'c)','Units','normalized','FontWeight','Bold','FontSize',11)