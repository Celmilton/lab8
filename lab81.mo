model lab81
Real p_cr = 41; 
Real N = 90;
Real q = 1; 
Real tau1 = 29; 
Real tau2 = 26;
Real p1 = 12.5;
Real p2 = 10.5; 
Real M1;
Real M2;

Real a1 = p_cr/(tau1*tau1*p1*p1*N*q);
Real a2 = p_cr/(tau2*tau2*p2*p2*N*q);
Real b = p_cr/(tau1*tau1*tau2*tau2*p1*p1*p2*p2*N*q);
Real c1 = (p_cr-p1)/(tau1*p1);
Real c2 = (p_cr-p2)/(tau2*p2);

initial equation
M1 = 7.4;
M2 = 8.4;

equation
der (M1) = (c1/c1)*M1-(a1/c1)*M1*M1-(b/c1)*M1*M2;
der(M2) = (c2/c1)*M2-(a2/c1)*M2*M2-(b/c1)*M1*M2;
annotation(experiment( StartTime = 0, StopTime = 30, Tolerance = 1e-06, interval = 0.05));

end lab81;
