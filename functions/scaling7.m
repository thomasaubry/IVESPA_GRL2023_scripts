function H=scaling7(parameters,inputs)
%function for the Aubry et al. (2017) scaling
%H= height (km a.v.l.)
%M=Mass eruption rate (kg/s)
%N=Brunt Vaisala (s^-1)
%W= wind regime parameter for the Aubry et al. (2017) scaling 
%pi = model parameter
M=inputs(:,1);
N=inputs(:,2);
W=inputs(:,3);
p1=parameters(1);
p2=parameters(2);

H=p1*(M.^0.25).*(N.^(-0.75))./(1+p2*W).^0.5;


end