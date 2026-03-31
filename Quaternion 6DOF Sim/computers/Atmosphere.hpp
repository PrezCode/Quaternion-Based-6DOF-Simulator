//Standard Atmosphere calculations from "Missile Guidance and Control Systems" by George M. Siouris, pg. 605
//Complex Atmosphere calculations from https://walter.bislins.ch/bloge/index.asp?page=Deriving+Equations+for+Atmospheric+Pressure+and+Density
#pragma once
#include <iostream>
#include <cmath>
#include <string>

class Atmosphere{
    public:
    Atmosphere(){}
    void generateStandardAtmosphere(double altitude){//Only good up to ~300,000ft
        if(altitude <= 11000.0){  
            temperatureAtmospheric = temperatureASL - 0.006499708*altitude;
            pressureAtmospheric = pressureASL*pow(1 - (2.255692257E-5)*altitude, 5.2561);
        }else if(altitude > 11000.0 && altitude <= 25000.0){
            temperatureAtmospheric = 216.66666667;
            pressureAtmospheric = pressureASL*(0.223358)*exp((-1.576883202E-4)*(altitude - 11000.0));
        }else if(altitude > 25000.0){
            temperatureAtmospheric = 216.66666667 + (3.000145816E-3)*(altitude - 25000);
            pressureAtmospheric = 2489.773467*exp((-1.576883202E-4)*(altitude - 25000.0));
        }
        airDensity = pressureAtmospheric/(R*temperatureAtmospheric);
    }
    void generateComplexAtmosphere(double altitude){//Should be good above 300,000ft
        selectAtmosphericLayer(-altitude);
        //Isothermic Layers
        if(Tgradient == 0.0){
            pressureAtmospheric = Pref*exp(-g/(R*Tref)*(-altitude - altitudeRef));
            temperatureAtmospheric = Tref;
            airDensity = densityRef*exp(-g/(R*Tref)*(-altitude - altitudeRef));
        }else if(Tgradient != 0.0 && altitude > -700000.0){//Gradient Layers
            double beta = g/(R*Tgradient);
            temperatureAtmospheric = Tgradient*(-altitude - altitudeRef) + Tref;
            pressureAtmospheric = Pref*pow(1 + (Tgradient*(-altitude - altitudeRef))/Tref, -beta);
            airDensity = densityRef*pow(1 + Tgradient*(-altitude - altitudeRef)/Tref, -beta - 1);
        }else if(Tgradient == 0.0 && altitude > -700000.0){//Exo-Atmosphere (SPACE)
            temperatureAtmospheric = 2676.95;
            pressureAtmospheric = 0;
            airDensity = 0.0;
        }
    }
    
    double getAirPressure(){return pressureAtmospheric;}
    double getAirTemperature(){return temperatureAtmospheric;}
    double getAirDensity(){return airDensity;}
    double getSoundBarrier(){return 20.037673*sqrt(temperatureAtmospheric);}
    private:
    void selectAtmosphericLayer(double altitude){
        if(altitude < 11000.0){//Block 0
            Tref = 288.15;
            Pref = 101325.0;
            Tgradient = -0.0065;
            altitudeRef = 0.0;
            densityRef = 1.225;
        }else if(altitude >= 11000.0 && altitude < 20000.0){//Block 1
            Tref = 216.65;
            Pref = 22632.1;
            Tgradient = 0.0;
            altitudeRef = 11000.0;
            densityRef = 0.36392;
        }else if(altitude >= 20000.0 && altitude < 32000.0){//Block 2
            Tref = 216.65;
            Pref = 5474.89;
            Tgradient = 0.001;
            altitudeRef = 20000.0;
            densityRef = 0.088035;
        }else if(altitude >= 32000.0 && altitude < 47000.0){//Block 3
            Tref = 228.65;
            Pref = 868.019;
            Tgradient = 0.0028;
            altitudeRef = 32000.0;
            densityRef = 0.013225;
        }else if(altitude >= 47000.0 && altitude < 51000.0){//Block 4
            Tref = 270.65;
            Pref = 110.9063;
            Tgradient = 0.0;
            altitudeRef = 47000.0;
            densityRef = 0.00142753;
        }else if(altitude >= 51000.0 && altitude < 71000.0){//Block 5
            Tref = 270.65;
            Pref = 66.9389;
            Tgradient = -0.0028;
            altitudeRef = 51000.0;
            densityRef = 0.0008616;
        }else if(altitude >= 71000.0 && altitude < 84852.0){//Block 6
            Tref = 214.65;
            Pref = 3.95642;
            Tgradient = -0.002;
            altitudeRef = 71000.0;  
            densityRef = 6.4211E-5;
        }else if(altitude >= 84852.0 && altitude < 100000.0){//Block 7
            Tref = 186.946;
            Pref = 0.37339868;
            Tgradient = 0.0;
            altitudeRef = 84852.0; 
            densityRef = 1.86347E-5;
        }else if(altitude >= 100000.0 && altitude < 110000.0){//Beginning of Exo-atmospheric Layers
            Tref = 186.946;
            Pref = 3.0075E-2;
            Tgradient = 0.005;
            altitudeRef = 100000.0;
            densityRef = 4.3681E-7;
        }else if(altitude >= 110000.0 && altitude < 120000.0){
            Tref = 236.946;
            Pref = 7.3544E-3;
            Tgradient = 0.01;
            altitudeRef = 110000.0;
            densityRef = Pref/(R*Tref);
        }else if(altitude >= 120000.0 && altitude < 150000.0){
            Tref = 336.946;
            Pref = 2.5217E-3;
            Tgradient = 0.02;
            altitudeRef = 120000.0;
            densityRef = Pref/(R*Tref);
        }else if(altitude >= 150000.0 && altitude < 160000.0){
            Tref = 936.946;
            Pref = 5.0617E-4;
            Tgradient = 0.015;
            altitudeRef = 150000.0;
            densityRef = Pref/(R*Tref);
        }else if(altitude >= 160000.0 && altitude < 170000.0){
            Tref = 1086.95;
            Pref = 3.6943E-4;
            Tgradient = 0.01;
            altitudeRef = 160000.0;
            densityRef = Pref/(R*Tref);
        }else if(altitude >= 170000.0 && altitude < 190000.0){
            Tref = 1186.95;
            Pref = 7.3544E-3;
            Tgradient = 0.007;
            altitudeRef = 170000.0;
            densityRef = Pref/(R*Tref);
        }else if(altitude >= 190000.0 && altitude < 230000.0){
            Tref = 1326.95;
            Pref = 1.6852E-4;
            Tgradient = 0.005;
            altitudeRef = 190000.0;
            densityRef = Pref/(R*Tref);
        }else if(altitude >= 230000.0 && altitude < 300000.0){
            Tref = 1526.95;
            Pref = 6.9604E-5;
            Tgradient = 0.004;
            altitudeRef = 230000.0;
            densityRef = Pref/(R*Tref);
        }else if(altitude >= 300000.0 && altitude < 400000.0){
            Tref = 1806.95;
            Pref = 1.8838E-5;
            Tgradient = 0.0033;
            altitudeRef = 300000.0;
            densityRef = Pref/(R*Tref);
        }else if(altitude >= 400000.0 && altitude < 500000.0){
            Tref = 2136.95;
            Pref = 4.0304E-6;
            Tgradient = 0.0026;
            altitudeRef = 400000.0;
            densityRef = Pref/(R*Tref);
        }else if(altitude >= 500000.0 && altitude < 600000.0){
            Tref = 2396.95;
            Pref = 1.0957E-6;
            Tgradient = 0.0017;
            altitudeRef = 500000.0;
            densityRef = Pref/(R*Tref);
        }else if(altitude >= 600000.0 && altitude < 700000.0){
            Tref = 2566.95;
            Pref = 3.4502E-7;
            Tgradient = 0.0011;
            altitudeRef = 600000.0;
            densityRef = Pref/(R*Tref);
        }else if(altitude > 700000.0){
            Tref = 2676.95;
            Pref = 0.0;
            Tgradient = 0.0;
            altitudeRef = 700000.0;
            densityRef = 0.0;
        }
    }
    double m2ft{3.28084}, ft2m{0.3048},
    altitude{0.0/*meters*/}, altitudeRef{0.0}, 
    airDensity{0.0/*kg/m^3*/}, densityRef{0.0},
    temperatureASL{288.1667/*Kelvin*/}, temperatureAtmospheric{0.0}, Tref{0.0/*K*/}, Tgradient{0.0/*K/km*/},
    pressureASL{101314.628/*Pa*/}, pressureAtmospheric{0.0}, Pref{0.0/*Pa*/},
    R{287.058/*gas constant*/}, 
    k{1.4/*specific heat for air*/},
    g{9.80665/*m/s^2*/};
};