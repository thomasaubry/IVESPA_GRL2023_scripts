function H=scaling4(parameters,inputs)
%function for the Hewett et al. (1971) scaling
%H = height (km a.v.l.)
%M =Mass eruption rate (kg/s)
%N =Brunt Vaisala (s^-1)
%W = Wind speed (m/s)
%pi = model parameter
p1=parameters(1);
M=inputs(:,1);
N=inputs(:,2);
W=inputs(:,3);

H=p1*(M.^(1/3)).*(N.^(-2/3)).*(W.^(-1/3));


end