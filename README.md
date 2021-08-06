# The SIR Model for Spread of Disease and its analysis for COVID-19

## Overview
Final project for my Bachelor's degree in Mathematics for Engineering at Politecnico di Torino. </br>
The aim of the thesis is to show how it is possible to model the spread of COVID-19 in Italy using a SIRD model. </br>
Before looking at the code I recommend checking out the available pdf *SIR_Theory_and_Simulations*.

## Part 1: Compartmental Models 

### The SIR Model

There are various approaches to modelling infectious diseases. The SIR model is a compartmental model, that requires the population to be divided into 3 categories: </br>
- Suscepitbles: individuals that can contract the disease;
- Infected: individuals who have contracted the disease;
- Removed: individuals that have recovered or died from the disease;

We also define two rates:  
- $\alpha$, rate of infection;
-  $\gamma$, rate of removal;

This allows us to express the dynamics of the spread of the disease by means of a differential equation model (please see the pdf previously mentioned for more insight on this).

This model was used to model an influenza outbreak that took place in an English boarding School in 1978, and gives very appreciable fitting.</br> T

### The SIRD Model

Since the removed category does not allow to differentiate between those who have recovered from the disease and those who have dies because of it, we can extend the SIR model to include such distinction. Specifically, we define:
- Recovered: infidividuals who have recovered from the disease;
- Dead: individuals who have died from the disease;

Moreover, we have:
- $\gamma$, recovery rate;
- $\delta$, death rate;

In simple cases (such as the influenza outbreak previously addressed), these rates can be considered as constants. However in more complicated cases, such as Covid-19, we make these parameters time-dependant.

## Part 2: Modelling the spread of Covid-19 in Italy

To model the spread of Covid-19 in Italy the SIRD model was used and data up to July 2020 was considered. For this reason the recovery rate $\gamma$ was considered constant, since there were no vaccines availbale at the time. As far as the infection rate is concerned instead, it was considered to be constant up to when draconian measures were imposed. From then on it is modeled as an exponentially decreasing function, since lockdown measures lead to a decrease in the transmission rate. The death rate was also modeled in this manner, since at the begining the mortality rate is higher given that the more severe cases are detected, while further on also milder cases are added. </br>

After having chosen the model to use, we proceed to solve an optimisation problem in order to find the rates that give the best possible fits. To address this problem we made use of the *fmincon* function provided by *MATLAB*, which is gradient based.

Subsequently, to highlight how this model is suited only for modelling the spread of the virus given all data, but not for making prediction, we apply the method considering only the data of the first 40 days. The performance decreses significantly (once again, see the pdf for details). We can imporve it slightly by means of the Tukey Loss function instead of an L2 Loss (the residuals for computing such losses are given by the difference between the real value for a given category and the predicted value obtained by using the SIRD model).

## Files in the repository

- *SIR_Theory_and_Simulations.pdf*: contains a detailed explanation of the thesis and the implemented mathematical methods;
- *plotSIR.m*: code to plot the curves for a SIR model given the initial state of the system and the various rates (which are supposed to be constant);
- *OptimisationExplained.m*: allows the user to better understand the use of the fmincon algorithm to obtain the optimal values for the rates by means of an example. Plots are provided to visually see the difference between the model with initial choice of rates vs the optimal choice;
- *COVID_SIRD_dataFitting.m*: models the spread of Covid-19 in Italy using the SIRD model and fmincon. In this case the rates are time-dependant (see Part 2);
- *COVID_SIRD_prediction.m*: prediction of the spread of Covid-19 in Italy using the SIRD model and fmincon. Here also the rates are time-dependant;
- *COVID_SIRD_prediction_Tukey.m*: like above, with the use of he Tukey Loss function instead of an L2 loss.
- *DataSIRD_COVID-19_Italy.xlsx* and *DataSIRExample.xlsx*: contain the data used to test the models.
