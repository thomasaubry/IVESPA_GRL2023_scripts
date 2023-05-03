%This script is just opening, reading and preprocessing some IVESPA data

%get the text and numerical data from the excel spreadsheet
[data textdata]=xlsread('IVESPA_GRL2023_data.xlsx');

%get the headers
headerln1=textdata(1,:);
headerln2=textdata(2,:);


%everything later follows a similar structure as the two lines below: first
%find the row corresponding to a target variable, then go read it and store
%it in a variable
i=find(strcmp(headerln2,'GVP volcano number'))-1;
IV_volcano=data(:,i);


i=find(strcmp(headerln2,'GVP eruption number'))-1;
IV_eruption=data(:,i);

i=find(strcmp(headerln2,'Vent altitude (m a.s.l.)'))-1;
IV_ventalt=data(:,i)/1000;%Here note I convert the vent height to km
i=find(strcmp(headerln2,'Eruption style'));
IV_style=textdata(3:end,i);
i=find(strcmp(headerln2,'IVESPA ID'));
IV_ivid=textdata(3:end,i);
i=find(strcmp(headerln2,'Plume morphology'));
IV_morpho=textdata(3:end,i);
i=find(strcmp(headerln2,'Ash Plume Top Method'));
IV_heighttech=textdata(3:end,i);

%Here note I convert all the plume heights from km above sea level to km
%above vent level
i=find(strcmp(headerln2,'Ash Plume Top (km asl) Best estimate'))-1;
IV_Htop_BE=data(:,i)-IV_ventalt;
i=find(strcmp(headerln2,'Ash Plume Top Best Estimate Flag'))-1;
IV_Htop_BE_flag=data(:,i);
i=find(strcmp(headerln2,'Ash Plume Top (km asl) Uncertainty'))-1;
IV_Htop_U=data(:,i);
i=find(strcmp(headerln2,'Ash Plume Top Uncertainty Flag'))-1;
IV_Htop_U_flag=data(:,i);
i=find(strcmp(headerln2,'Spreading height of ash (km asl) Best estimate'))-1;
IV_Hspread_BE=data(:,i)-IV_ventalt;
i=find(strcmp(headerln2,'Spreading height of ash Best Estimate Flag'))-1;
IV_Hspread_BE_flag=data(:,i);
i=find(strcmp(headerln2,'Spreading height of ash (km asl) Uncertainty'))-1;
IV_Hspread_U=data(:,i);
i=find(strcmp(headerln2,'Spreading height of ash Uncertainty Flag'))-1;
IV_Hspread_U_flag=data(:,i);
i=find(strcmp(headerln2,'SO2 height (km asl) Best estimate'))-1;
IV_Hso2_BE=data(:,i)-IV_ventalt;
i=find(strcmp(headerln2,'SO2 height Best Estimate Flag'))-1;
IV_Hso2_BE_flag=data(:,i);
i=find(strcmp(headerln2,'SO2 height (km asl) Uncertainty'))-1;
IV_Hso2_U=data(:,i);
i=find(strcmp(headerln2,'SO2 height Uncertainty Flag'))-1;
IV_Hso2_U_flag=data(:,i);
i=find(strcmp(headerln2,'Mass of PDC (kg) Best estimate'))-1;
IV_PDC_BE=data(:,i);
i=find(strcmp(headerln2,'Mass of PDC Uncertainty'))-1;
IV_PDC_U=data(:,i);
i=find(strcmp(headerln2,'Exit velocity (m/s) Best estimate'))-1;
IV_u0_BE=data(:,i);
i=find(strcmp(headerln2,'Exit velocity (m/s) Uncertainty'))-1;
IV_u0_U=data(:,i);


i=find(strcmp(headerln2,'Magma water content (wt.%) Best estimate'))-1;
IV_n0_BE=data(:,i)/100;%I convert the mass fraction from wt.% to fraction
i=find(strcmp(headerln2,'Magma water content (wt.%) Uncertainty'))-1;
IV_n0_U=data(:,i)/100;

i=find(strcmp(headerln2,'Magma temperature (deg C) Best estimate'))-1;
IV_T0_BE=data(:,i)+273.15;%I convert the temperature from C to K
i=find(strcmp(headerln2,'Magma temperature (deg C) Uncertainty'))-1;
IV_T0_U=data(:,i);


%==========================================================================
%below are various atmospheric parameters from two family of reanalyses

i=find(strcmp(headerln2,'Average wind speed (m/s, NOAA reanalyses)'))-1;
IV_W_NOAA=data(:,i);
i=find(strcmp(headerln2,'Average Brunt-Väisälä frequency (/s, NOAA reanalyses)'))-1;
IV_N_NOAA=data(:,i);
i=find(strcmp(headerln2,'Average wind speed (m/s, ERA reanalyses)'))-1;
IV_W_ERA=data(:,i);
i=find(strcmp(headerln2,'Average Brunt-Väisälä frequency (/s, ERA reanalyses)'))-1;
IV_N_ERA=data(:,i);
i=find(strcmp(headerln2,'Average wind shear (s-1, ERA reanalyses)'))-1;
IV_Wshear_ERA=data(:,i);
i=find(strcmp(headerln2,'Average wind shear (s-1, NOAA reanalyses)'))-1;
IV_Wshear_NOAA=data(:,i);
i=find(strcmp(headerln2,'Average tropopause height (km, ERA reanalyses)'))-1;
IV_tph_ERA=data(:,i);
i=find(strcmp(headerln2,'Average tropopause height (km, NOAA reanalyses)'))-1;
IV_tph_NOAA=data(:,i);

i=find(strcmp(headerln2,'Average relative humidity (%, NOAA reanalyses)'))-1;
IV_RH_NOAA=data(:,i);
i=find(strcmp(headerln2,'Average relative humidity (%, ERA reanalyses)'))-1;
IV_RH_ERA=data(:,i);


i=find(strcmp(headerln2,'Relative humidity at vent altitude (%, NOAA reanalyses)'))-1;
IV_RH0_NOAA=data(:,i);
i=find(strcmp(headerln2,'Relative humidity at vent altitude (%, ERA reanalyses)'))-1;
IV_RH0_ERA=data(:,i);
i=find(strcmp(headerln2,'Temperature at vent altitude (K, NOAA reanalyses)'))-1;
IV_T0_NOAA=data(:,i);
i=find(strcmp(headerln2,'Temperature at vent altitude (K, ERA reanalyses)'))-1;
IV_T0_ERA=data(:,i);
i=find(strcmp(headerln2,'Pressure at vent altitude (Pa, NOAA reanalyses)'))-1;
IV_P0_NOAA=data(:,i);
i=find(strcmp(headerln2,'Pressure at vent altitude (Pa, , ERA reanalyses)'))-1;
IV_P0_ERA=data(:,i);


%below are the values that atmospheric parameters would have for each event
%if that event occured with the same vent and plume altitude, but with
%atmospheric profiles equal to the average of all atmospheric profiles in
%IVESPA

i=find(strcmp(headerln2,'Average wind speed (m/s, reference profile)'))-1;
IV_W_ref=data(:,i);
i=find(strcmp(headerln2,'Average Brunt-Väisälä frequency (/s, reference profile)'))-1;
IV_N_ref=data(:,i);
i=find(strcmp(headerln2,'Average wind shear (s-1, reference profile)'))-1;
IV_Wshear_ref=data(:,i);
i=find(strcmp(headerln2,'Average relative humidity (%, reference profile)'))-1;
IV_RH_ref=data(:,i);


%For atmospheric parameters, I define the best estimate as the average of
%the ERA and NOAA values
IV_N_BE=0.5*(IV_N_ERA+IV_N_NOAA);
IV_W_BE=0.5*(IV_W_ERA+IV_W_NOAA);
IV_N_U=0.5*abs(IV_N_ERA-IV_N_NOAA);
IV_W_U=0.5*abs(IV_W_ERA-IV_W_NOAA);
IV_Wshear_BE=0.5*(IV_Wshear_ERA+IV_Wshear_NOAA);
IV_RH_BE=0.5*(IV_RH_ERA+IV_RH_NOAA);
IV_RHA0_BE=0.5*(IV_RH0_ERA+IV_RH0_NOAA);
IV_TA0_BE=0.5*(IV_T0_ERA+IV_T0_NOAA);
IV_PA0_BE=0.5*(IV_P0_ERA+IV_P0_NOAA);

%==========================================================================


i=find(strcmp(headerln2,'Median grain size (phi unit)'))-1;
IV_medianphi_BE=data(:,i);
i=find(strcmp(headerln2,'Latitude (deg. N)'))-1;
IV_latitude=data(:,i);

% Here note I convert the duration to second
i=find(strcmp(headerln2,'Duration Best estimate (hours)'))-1;
IV_duration_BE=data(:,i)*3600;
i=find(strcmp(headerln2,'Duration Best Estimate Flag'))-1;
IV_duration_BE_flag=data(:,i);
i=find(strcmp(headerln2,'Duration Uncertainty (hours)'))-1;
IV_duration_U=data(:,i)*3600;
i=find(strcmp(headerln2,'Duration Uncertainty Flag'))-1;
IV_duration_U_flag=data(:,i);

i=find(strcmp(headerln2,'TEM Best estimate (kg)'))-1;
IV_TEM_BE=data(:,i);
i=find(strcmp(headerln2,'TEM Best Estimate Flag'))-1;
IV_TEM_BE_flag=data(:,i);
i=find(strcmp(headerln2,'TEM Uncertainty (low; kg)'))-1;
IV_TEM_UL=data(:,i);
i=find(strcmp(headerln2,'TEM Uncertainty (high; kg)'))-1;
IV_TEM_UU=data(:,i);
i=find(strcmp(headerln2,'TEM Uncertainty Flag'))-1;
IV_TEM_U_flag=data(:,i);

%==========================================================================
%Last, below I calculate some regime parameters for scaling that include a
%dependence on wind. Details are not provided in this script, please refer
%to the supporting information of our GRL 2023 paper.


%Aubry et al. 2017
n0_filled=IV_n0_BE;n0_filled(isnan(n0_filled))=mean(n0_filled(~isnan(n0_filled)));
T0_filled=IV_T0_BE;T0_filled(isnan(T0_filled))=mean(T0_filled(~isnan(T0_filled)));
IV_Ws_BE=IV_W_BE./(1.85*(461.5*n0_filled.*T0_filled).^0.5);

%Woodhouse et al 2013
IV_MER_BE=IV_TEM_BE./IV_duration_BE;
IV_Riw_BE=IV_Wshear_BE./IV_N_BE;

%Degruyter and Bonadonna 2012
IV_Vs_BE=IV_W_BE./(IV_N_BE.*IV_MER_BE).^0.25;
IV_VsDeg_BE=IV_W_BE./(1000*IV_Htop_BE.*IV_N_BE);
IV_VsDeg_BE(isnan(IV_VsDeg_BE))=IV_W_BE(isnan(IV_VsDeg_BE))./(1000*IV_Hspread_BE(isnan(IV_VsDeg_BE)).*IV_N_BE(isnan(IV_VsDeg_BE)));
IV_VsDeg_BE(isnan(IV_VsDeg_BE))=IV_W_BE(isnan(IV_VsDeg_BE))./(1000*IV_Hso2_BE(isnan(IV_VsDeg_BE)).*IV_N_BE(isnan(IV_VsDeg_BE)));


%To calculate their PI parameter, I assume beta=0.5 and alpha = 0.1
IV_PI_BE=(1000*IV_Htop_BE.*IV_N_BE)./IV_W_BE*(0.1^2/(1.8*0.5^2));
IV_PI_BE(isnan(IV_PI_BE))=(1000*IV_Hso2_BE(isnan(IV_PI_BE)).*IV_N_BE(isnan(IV_PI_BE)))./IV_W_BE(isnan(IV_PI_BE))*(0.1^2/(1.8*0.5^2));
IV_PI_BE(isnan(IV_PI_BE))=(1000*IV_Hspread_BE(isnan(IV_PI_BE)).*IV_N_BE(isnan(IV_PI_BE)))./IV_W_BE(isnan(IV_PI_BE))*(0.1^2/(1.8*0.5^2));
IV_PI_U=IV_PI_BE.*((0.5*IV_Htop_U./IV_Htop_BE).^2+(0.5*IV_N_U./IV_N_BE).^2+(0.5*IV_W_U./IV_W_BE).^2+2*(0.025./0.1).^2+2*(0.2./0.5).^2).^0.5;
IV_PI_U(isnan(IV_PI_U))=IV_PI_BE(isnan(IV_PI_U)).*((0.5*IV_Hso2_U(isnan(IV_PI_U))./IV_Hso2_BE(isnan(IV_PI_U))).^2+(0.5*IV_N_U(isnan(IV_PI_U))./IV_N_BE(isnan(IV_PI_U))).^2+(0.5*IV_W_U(isnan(IV_PI_U))./IV_W_BE(isnan(IV_PI_U))).^2+2*(0.025./0.1).^2+2*(0.2./0.5).^2).^0.5;
IV_PI_U(isnan(IV_PI_U))=IV_PI_BE(isnan(IV_PI_U)).*((0.5*IV_Hspread_U(isnan(IV_PI_U))./IV_Hspread_BE(isnan(IV_PI_U))).^2+(0.5*IV_N_U(isnan(IV_PI_U))./IV_N_BE(isnan(IV_PI_U))).^2+(0.5*IV_W_U(isnan(IV_PI_U))./IV_W_BE(isnan(IV_PI_U))).^2+2*(0.025./0.1).^2+2*(0.2./0.5).^2).^0.5;


%==========================================================================
%A few last random parameters and some cleaning
i=find(strcmp(headerln2,'Region or volcano'));
IV_region=textdata(3:end,i);
i=find(strcmp(headerln2,'Event Year'))-1;
IV_year=data(:,i);

%below I just delete variables I don't need
clear i data textdata headerln1 headerln2

%et voila!
