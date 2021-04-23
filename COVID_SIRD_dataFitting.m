% To understand the following code we recommend reading the
% "SIR_Theory_and_Simulations" pdf,
% specifically slides 17, 18, 19 and 20.

%% Get Data from Excel File
[Data, Days, raw]=xlsread('DataSIRD_COVID-19_Italy.xlsx');
S_exp=Data(:,1)';
I_exp=Data(:,2)';
R_exp=Data(:,5)';
D_exp=Data(:,6)';
N=Data(1,7);
date=datetime(Days(3:165,1),'Format','dd/MM/yyyy');

t=linspace(0,length(I_exp),length(I_exp)+1); %from 20th Feb to 31July
x0=ones(1,5); %x=[alpha0 taualpha delta0 taudelta gamma]
lb=[0.0,10,0.0,10,1/100]; %lower bound for x
ub=[0.55,20, 0.5,30,1/10]; %upper bound for x

I0=I_exp(1);
R0=R_exp(1);
S0 = N-I0-R0-D_exp(1);

%% Gradient based optimisation
options = optimoptions(@fmincon,'Display', 'off','PlotFcn', { @optimplotfval });
[x,fval] = fmincon(@(x)sird_covid_ols(x,t,S0,I0,R0,S_exp,I_exp,R_exp,D_exp,N),x0,[],[],[],[],lb,ub,[],options);

%% Solving the model with the best fit parameters
t=linspace(0,length(S_exp)-1,length(S_exp));
[t,I,R,D,S]=sird_covid_solver(x,t,S0,I0,R0,N);

%% Define parameters of the outbreak
alpha_t=ones(length(t));
delta_t=alpha_t;

for i=1:length(t)
    if i<=19 %Days before the lockdown
        alpha_t(i)=x(1);
        delta_t(i)=x(3);
    else
        alpha_t(i)=x(1)*exp((-t(i)+19)/x(2));
        delta_t(i)=x(3)*exp((-t(i)+19)/x(4));
        
    end
end

gamma=x(5);
RN=(alpha_t'./(gamma+delta_t'));
res=sird_covid_ols(x,t,S0,I0,R0,S_exp,I_exp,R_exp,D_exp,N);

%% Plot S
figure
set(gcf,'color','w');
plot(date,S_exp,'go',date,S([1:length(S_exp)]),'-g')
title('Susceptible')
xlabel('Days')
ylabel('Individuals')
ylim([5.97*10^7 6.02*10^7])
legend('Data S','Fit S','Location','east');
box on;

%% Plot I
figure
set(gcf,'color','w');
plot(date,I_exp,'ro',date,I([1:length(I_exp)]),'-r')
title('Infected')
xlabel('Days')
ylabel('Individuals')
ylim([0 17*10^4])
legend('Data I','Fit I','Location','east');

%% Plot R
figure
set(gcf,'color','w');
plot(date,R_exp,'bo',date,R([1:length(I_exp)]),'-b')
title('Recovered')
xlabel('Days')
ylabel('Individuals')
legend('Data R','Fit R','Location','east');

%% Plot D
figure
set(gcf,'color','w');
plot(date,D_exp,'ko',date,D([1:length(I_exp)]),'-k')
title('Deaths')
xlabel('Days')
ylabel('Individuals')
legend('Data D','Fit D','Location','east');

%% Plot Together
figure
set(gcf,'color','w');
data_long=dateshift(date(1),'start','day',0:t(end));
plot(date,I_exp,'ro',data_long,I,'-r',date,R_exp,'bo',data_long,R,'-b',date,D_exp,'ko',data_long,D,'-k')
ylabel('Individuals')
xlabel('Date (dd/MM)')
xlim(datetime(2020,[2 7],[15 31]))
ylim([0 2.2*10^5])
datetick('x','dd-mm','keeplimits')
legend('Data I','Fit I', 'Data R','Fit R', 'Data D', 'Fit D', 'Location', 'east')
title('SIRD Model for COVID-19 in Italy')
box on;

function [t,I,R,D,S]=sird_covid_solver(x,t,S0,I0,R0,N)

alpha0=x(1);
tau_alpha=x(2);
delta0=x(3);
tau_delta=x(4);
gamma0=x(5);

options=odeset('NonNegative',(1:3));
[t,y] = ode45(@(t,y)sird_covid_odes_ita(t,y,alpha0,tau_alpha,delta0,tau_delta,gamma0,N),t,[S0; I0; R0],options);

S=y(:,1);
I=y(:,2);
R=y(:,3);
D=N-S-I-R;
end

function ols_t=sird_covid_ols(x,t,S0,I0,R0,S_exp,I_exp,R_exp,D_exp,N)

[t,I,R,D,S]=sird_covid_solver(x,t,S0,I0,R0,N);

ols_s=norm(S_exp-S,2)^2;
ols_i=norm(I_exp-I,2)^2;
ols_r=norm(R_exp-R,2)^2;
ols_d=norm(D_exp-D,2)^2;

ols_t=ols_s+ols_i+ols_r+ols_d;
end

function dydt = sird_covid_odes_ita(t,y,alpha0,tau_alpha,delta0,tau_delta,gamma0,N)

 if t<=19 % Before the lockdown
     alpha=alpha0;
     delta=delta0;
 else
    alpha=alpha0*exp((-t+15)/(tau_alpha));
    delta=delta0*exp((-t+15)/tau_delta);

 end
gamma=gamma0;

S=y(1);
I=y(2);
R=y(3);

dS_dt=-alpha*(I/N)*S;
dI_dt=alpha*(I/N)*S-gamma*I-delta*I;
dR_dt=gamma*I;
dydt = [dS_dt; dI_dt; dR_dt];

end
