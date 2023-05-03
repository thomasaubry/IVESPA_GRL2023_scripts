function H=scaling1(parameters,inputs)
%function for the "canonical" power-law model
%H= height (km a.v.l.)
%M=Mass eruption rate (kg/s)
%pi = model parameter
M=inputs(:,1);
p1=parameters(1);
p2=parameters(2);
H=p1*M.^p2;


end