function H=scaling3(parameters,inputs)
%function for the Morton et al. (1956) scaling power-law model
%H= height (km a.v.l.)
%M=Mass eruption rate (kg/s)
%N=Brunt Vaisala (s^-1)
%pi = model parameter

p1=parameters(1);
M=inputs(:,1);
N=inputs(:,2);
H=p1*(M.^0.25).*(N.^(-0.75));


end