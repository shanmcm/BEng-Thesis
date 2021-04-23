function plotSIR(alpha, gamma, N, Io, period)

    %  This function gives a visual representation of the SIR model. Specifically, it plots the curves of
    %  Susceptibles, Infected and Removed individuals over the course of time.
    %  The input arguments are the rates of transmission (alpha) and removal (gamma), the total population (N), the initial number of infected individuals (Io) and the period of time considered in days (period).
    %  Here is an example of possble input parameters:
    %
    %         N=763;
    %         alpha=0.0026;
    %         gamma=0.5;
    %         Io=1;
    %         period=15;
    %
    %  Just run the following in the command window: plotSIR(0.0026, 0.5, 763, 1, 15)
   
    %Data
    t=linspace(0,period);
    fit = SIR([alpha gamma]', t);
    R0=N*alpha/gamma;
    
    %Plot
    figure
    set(gcf,'color','w');
    hold on
    plot(t,fit(:,1),'g');
    plot(t,fit(:,2),'r');
    plot(t,fit(:,3),'b');
    grid on;
    xlabel('Days'); 
    ylabel('Number of individuals');
    legend('S','I','R');
    title([' \alpha = ' num2str(alpha) ', \gamma = ' num2str(gamma) ', R0 = ' num2str(R0)]);
    hold off

    function D=SIR(parameters,t) %paramters contains beta and gamma
        
        zo=[N-Io, Io, 0]'; %Initial conditions
        [T,P]=ode15s(@DifEq,t,zo); 
        
            function dZ=DifEq(t,z)

                dzdt=zeros(3,1);
                %Equations of SIR model
                dzdt(1)=-parameters(1).*z(1).*z(2); %dS/dt
                dzdt(2)= parameters(1).*z(1).*z(2)-parameters(2).*z(2); %dI/dt
                dzdt(3)= parameters(2).*z(2); %dR/dt
                dZ=dzdt;

            end

        D=P;

    end

end
    