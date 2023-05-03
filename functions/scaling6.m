function H=scaling6(parameters,inputs)
%function for the Woodhouse et al. (2013) scaling
%H= height (km a.v.l.)
%M=Mass eruption rate (kg/s)
%N=Brunt Vaisala (s^-1)
%W= gradient Richardson number as defined in Woodhouse et al. (2013)
%we= intermediate fractional function of W (see Woodhouse et al. (2013)
%pi = model parameter
M=inputs(:,1);
N=inputs(:,2);
W=inputs(:,3);
p1=parameters(1);
p2=parameters(2);
we=(1+(0.87+0.50*p2)*W)./(1+(1.09+0.32*p2)*W+(0.06+0.03*p2)*W.^2);
H=p1*(M.^0.25).*(N.^(-0.75)).*we;


end