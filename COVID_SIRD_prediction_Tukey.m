% To understand the following code we recommend reading the
% "SIR_Theory_and_Simulations" pdf,
% specifically slides 17, 18, 19, 20, 21, 22 and 23.

%% Get Data from Excel File
[Data, Days, raw]=xlsread('DataSIRD_COVID-19_Italy.xlsx');
S_exp=Data(:,1)';
I_exp=Data(:,2)';
R_exp=Data(:,5)';
D_exp=Data(:,6)';
N=Data(1,7);
date=datetime(Days(3:165,1),'Format','dd/MM/yyyy');
x0=ones(1,5); %x=[alpha0 taualpha delta0 taudelta gamma]
lb=[0.0,10,0.0,10,1/100]; %lower bound for x
ub=[0.55,20, 0.5,30,1/10]; %upper bound for x
I0=I_exp(1);
R0=R_exp(1);
S0 = N-I0-R0-D_exp(1);

%% Define Data for Fitting
n=40;
S_fit=Data(1:n,1)';
I_fit=Data(1:n,2)';
R_fit=Data(1:n,5)';
D_fit=Data(1:n,6)';
t_fit=linspace(0,length(I_fit),length(I_fit)+1);

%% Gradient based optimisation
options = optimoptions(@fmincon,'Display', 'off','PlotFcn', { @optimplotfval });
[x,fval,exitflag,output,lambda,grad,hessian] = ...
fmincon(@(x)sird_covid_tukey(x,t_fit,S0,I0,R0,S_fit,I_fit,R_fit,D_fit,N),x0,[],[],[],[],lb,ub,[],options);

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

%% Reproduction number (it is not a constant, but a vector)
RN=(alpha_t'./(gamma+delta_t'));
prediction_end_date=datetime(['31/03/2020'; '31/03/2020'],'Format','dd/MM/yyyy');
res=sird_covid_tukey(x,t,S0,I0,R0,S_exp,I_exp,R_exp,D_exp,N);

%Plot S
figure
set(gcf,'color','w');
hold on
plot(date,S_exp,'go',date,S([1:length(S_exp)]),'-g')
prediction_end_value = [5.97*10^7; 6.01*10^7];
plot(prediction_end_date, prediction_end_value,'Color', [17 17 17]/255,'LineStyle','--')
title('Susceptible')
xlabel('Days')
ylabel('Individuals')
ylim(prediction_end_value)
legend('Data S','Fit S','Location','east');
box on;

%Plot I
figure
set(gcf,'color','w');
hold on
plot(date,I_exp,'ro',date,I([1:length(I_exp)]),'-r')
prediction_end_value = [0 1.3*10^5];
plot(prediction_end_date, prediction_end_value,'Color', [17 17 17]/255,'LineStyle','--')
title('Infected')
xlabel('Days')
ylabel('Individuals')
ylim(prediction_end_value)
legend('Data I','Fit I','Location','east');

%Plot R
figure
set(gcf,'color','w');
hold on
plot(date,R_exp,'bo',date,R([1:length(I_exp)]),'-b')
prediction_end_value = [0 2.2*10^5];
plot(prediction_end_date, prediction_end_value,'Color', [17 17 17]/255,'LineStyle','--')
title('Recovered')
xlabel('Days')
ylim(prediction_end_value)
ylabel('Individuals')
legend('Data R','Fit R','Location','east');

%Plot D
figure
set(gcf,'color','w');
hold on
plot(date,D_exp,'ko',date,D([1:length(I_exp)]),'-k')
prediction_end_value = [0 4.1*10^4];
plot(prediction_end_date, prediction_end_value,'Color', [17 17 17]/255,'LineStyle','--')
title('Deaths')
xlabel('Days')
ylim(prediction_end_value)
ylabel('Individuals')
legend('Data D','Fit D','Location','east');

%Plot Together
figure
set(gcf,'color','w');
hold on
data_long=dateshift(date(1),'start','day',0:t(end));
plot(date,I_exp,'ro',data_long,I,'-r',date,R_exp,'bo',data_long,R,'-b',date,D_exp,'ko',data_long,D,'-k')
prediction_end_value = [0; 2.2*10^5];
plot(prediction_end_date, prediction_end_value,'Color', [17 17 17]/255,'LineStyle','--')
ylabel('Individuals')
xlabel('Date (dd/MM)')
xlim(datetime(2020,[2 7],[15 31]))
ylim(prediction_end_value)
datetick('x','dd-mm','keeplimits')
legend('Data I','Fit I', 'Data R','Fit R', 'Data D', 'Fit D', 'Location', 'east')
title('SIRD Prediction for COVID-19 in Italy vs Data')
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

%% Robust Bisquare M-Estimator with Tukey function
function f=sird_covid_tukey(x,t,S0,I0,R0,S_fit,I_fit,R_fit,D_fit,N)

[t,I,R,D,S]=sird_covid_solver(x,t,S0,I0,R0,N);

rho_s=0;
rho_i=0;
rho_r=0;
rho_d=0;

cS=mad(S(1:length(S_fit))-S_fit')*4.685;
cI=mad(I(1:length(I_fit))-I_fit')*4.685; 
cR=mad(R(1:length(R_fit))-R_fit')*4.685;
cD=mad(D(1:length(D_fit))-D_fit')*4.685;

for i=1:length(S_fit)
    
    errs=S(i)-S_fit(i);
    if abs(errs)<=cS
        rho_s=rho_s+cS^2/6*(1-(1-(errs/cS)^2)^3);
    else
        rho_s=rho_s+cS^2/6;
    end
    
    erri=I(i)-I_fit(i);
    if abs(erri)<=cI
        rho_i=rho_i+cI^2/6*(1-(1-(erri/cI)^2)^3);
    else
        rho_i=rho_i+cI^2/6;
    end
    
    errr=R(i)-R_fit(i);
    if abs(errr)<=cR
        rho_r=rho_r+cR^2/6*(1-(1-(errr/cR)^2)^3);
    else
        rho_r=rho_r+cR^2/6;
    end
    
    errd=D(i)-D_fit(i);
    if abs(errd)<=cD
        rho_d=rho_d+cD^2/6*(1-(1-(errd/cD)^2)^3);
    else
        rho_d=rho_d+cD^2/6;
    end
    
end

f=rho_s+rho_i+rho_r+rho_d;
end

function dydt = sird_covid_odes_ita(t,y,alpha0,tau_alpha,delta0,tau_delta,gamma0,N)

 if t<=19 %Before lockdown
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
