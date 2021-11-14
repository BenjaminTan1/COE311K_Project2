%parameters:
%T=50; 
%g_vec and f_vec solved in main

x=delta;                                                 %vector of delta we want to optimize
y=tau;                                                   %vector of tau we want to optimize
x0=[0.001, 0.001, 0.001, 0.001];                         %initial guess for delta/drug volume fraction
y0=[9.7955, 19.5171, 30.6384, 40.6257];                  %initial guess for tal/time of ith drug input

xmin=[0, 0, 0, 0]; xmax=[0.01, 0.01, 0.01, 0.01];        %min and max values of delta

Ntreat=4;
ymin=zeros(1, 4);
ymax=zeros(1, 4);

for i=1:Ntreat                                          %for loop to creat tau min and max
    ymin(i)=5+10*(i-1);
    ymax(i)=15+10*(i-1);
end
 
tolx=1.e-9; tolfun=1.e-9; maxiter=400;                  %given

[xopt, yopt, Jval]=fmincon(J,x0,y0,xmin,xmax,ymin,ymax);
disp xopt
disp yopt
disp Jval
