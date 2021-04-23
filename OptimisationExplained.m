% To understand the following code we recommend reading the
% "SIR_Theory_and_Simulations" pdf,
% specifically slides 12 and 13.



%% Get Data from Excel File
Data=xlsread('DataSIRExample.xlsx');
S_exp=Data(:,1)';
I_exp=Data(:,2)';
R_exp=Data(:,3)';
N=Data(1,4);

t_exp=linspace(0,length(I_exp)-1,length(I_exp)); %from 20th Feb to 31July
x0=[0.0015 0.35]; %x=[alpha gamma], this is our initial guess of the values
lb=[0.0,1/100]; %lower bound for x
ub=[0.02,0.6]; %upper bound for x 

S0=S_exp(1);
I0=I_exp(1);
R0=R_exp(1);

S_fit=S_exp(1:4);
I_fit=I_exp(1:4);
R_fit=R_exp(1:4);
t_fit=linspace(0,length(I_fit)-1,length(I_fit));

%% Plot Initial Guess
t=linspace(0,15);
[t,S,I,R]=sird_covid_solver(x0,t,S0,I0,R0);

figure
set(gcf,'color','w');
plot(t_exp,S_exp,'go',t,S,'-g',t_exp,I_exp,'ro',t,I,'-r',t_exp,R_exp,'bo',t,R,'-b')
xlabel('Days')
ylabel('Individuals')
legend('Data S','Initial Fit S','Data I','Initial Fit I', 'Data R','Initial Fit R', 'Location', 'east')
title([' \alpha = ' num2str(x0(1)) ', \gamma = ' num2str(x0(2)) ', R0 = ' num2str(N*x0(1)/x0(2)) ]);

%% Initial value of function to be minimised
res_ini=sird_covid_ols(x0,t_exp,S0,I0,R0,S_exp,I_exp,R_exp);

%% Gradient based optimisation
options = optimoptions(@fmincon,'Display', 'off','PlotFcn', { @optimplotfval });
[x,fval] = fmincon(@(x)sird_covid_ols(x,t_exp,S0,I0,R0,S_exp,I_exp,R_exp),x0,[],[],[],[],lb,ub,[],options);

%% Solving the model with the best fit parameters
t=linspace(0,15);
[t,S,I,R]=sird_covid_solver(x,t,S0,I0,R0);

%% Define parameters of the outbreak
alpha=x(1);
gamma=x(2);
RN=N*alpha/gamma;

%% Plot Optimised
figure
hold on
set(gcf,'color','w');
plot(t_exp,S_exp,'go',t,S,'-g',t_exp,I_exp,'ro',t,I,'-r',t_exp,R_exp,'bo',t,R,'-b')
%plot([3 3],[0 500],'--k')
xlabel('Days')
ylabel('Individuals')
legend('Data S','Best Fit S','Data I','Best Fit I', 'Data R','Best Fit R', 'Location', 'east')
title([' \alpha = ' num2str(alpha) ', \gamma = ' num2str(gamma) ', R0 = ' num2str(RN) ]);

%% Value of minimised function
res_best=sird_covid_ols(x,t_exp,S0,I0,R0,S_exp,I_exp,R_exp);

function [t,S,I,R]=sird_covid_solver(x,t,S0,I0,R0)

alpha=x(1);
gamma=x(2);
[t,y] = ode15s(@(t,y)sird_covid_odes_ita(t,y,alpha,gamma),t,[S0; I0; R0]);

S=y(:,1);
I=y(:,2);
R=y(:,3);
end

function ols_t=sird_covid_ols(x,t,S0,I0,R0,S_exp,I_exp,R_exp)

[t,S,I,R]=sird_covid_solver(x,t,S0,I0,R0);

ols_s=norm(S_exp-S,2)^2;
ols_i=norm(I_exp-I,2)^2;
ols_r=norm(R_exp-R,2)^2;
ols_t=ols_s+ols_i+ols_r;

end

function dydt = sird_covid_odes_ita(t,y,alpha,gamma)

S=y(1);
I=y(2);
R=y(3);

dS_dt=-alpha*I*S;
dI_dt=alpha*I*S-gamma*I;
dR_dt=gamma*I;
dydt = [dS_dt; dI_dt; dR_dt];

end
