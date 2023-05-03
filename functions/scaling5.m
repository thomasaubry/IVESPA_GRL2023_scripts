function H=scaling5(parameters,inputs)
%function for the Degruyter and Bonadonna (2012) scaling
%H= height (km a.v.l.)
%M=Mass eruption rate (kg/s)
%N=Brunt Vaisala (s^-1)
%W= wind regime parameter for the Degruyter and Bonadonna (2012) scaling as
%defined in Aubry et al. (2017)
%we= intermediate fractional function of W (see Aubry et al. (2017))
%pi = model parameter
M=inputs(:,1);
N=inputs(:,2);
W=inputs(:,3);
p1=parameters(1);
p2=parameters(2);
we=(1+0.17*p2*W+0.00061*(p2*W).^2)./(1+0.48*p2*W+0.0072*(p2*W).^2);
H=p1*(M.^0.25).*(N.^(-0.75)).*we;


end