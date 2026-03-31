//Main file to test quaternion based simulator.
//Assume all values should be converted to metric and radians
#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include <iomanip>
#include <string>
#include "computers/QuaternionSimulation.hpp"
#include "models/NASA_ISS_B17.hpp"

void printInitialState(double* state){
        std::cout << std::setw(10) << " " 
        << std::setw(10) << "Xe" << std::setw(10) << "Ye" << std::setw(10) << "Ze" << std::setw(10) << "Ue" << std::setw(10) << "Ve" << std::setw(10) << "We" 
        << std::setw(10) << "P" << std::setw(10) << "Q" << std::setw(10) << "R"<< std::setw(10) << "Phi" << std::setw(10) << "Theta" << std::setw(10) << "Psi" 
        << std::endl;
        std::cout << std::setw(10) << "Initial: ";
        for(int i = 0; i < 12; ++i){std::cout << std::setw(10) << state[i] << " ";}
        std::cout << std::endl;
}

int main(){
    std::ofstream dataOut("MATLAB/Object_Simulation_Data.txt");
    std::array<double, 27> modelData{0.0};
    std::array<double, 12> finalState{0.0};
    std::string initial, final;
    Atmosphere air;
    double dt{0.01}, r2d{180.0/M_PI}, d2r{M_PI/180.0}, m2ft{3.28084}, ft2m{0.3048},
    initialState[12]{-4315967.74, 960356.20, 5167269.53, 129.091037, -7491.513855, 1452.515654, 0*d2r, -0.065*d2r, 0*d2r, 0*d2r, -11.60*d2r, 0*d2r};//Xe, Ye, Ze, Ue, Ve, We, P, Q, R, Phi, Theta, Psi; if spherical Earth, use decimal Lat-Long for X and Y
    int simEnd{6300}, simTimer{0};
    //Create and upload model
    InternationalSpaceStation ISS;
    ISS.getModel(modelData);
    QuaternionSimulator testObject(modelData, ObjectType::Miscellaneous, Earth::Spherical, Gravity::Calculated);
    testObject.setState(initialState);
    //Display initial state
    printInitialState(initialState);
    //Run simulation
    testObject.initialization();
    while(simTimer < simEnd/dt){
        testObject.iterate(dt);
        testObject.getState(finalState);
        dataOut << simTimer*dt << " ";
        for(int i = 0; i < 12; ++i){dataOut << finalState[i] << " ";}
        dataOut << std::endl;
        if(abs(finalState[2]) < 0.0){break;}
        simTimer++;
    }
    //Print final state
    std::cout << std::setw(10) << "Final: ";
    for(int value : finalState){std::cout << std::setw(10) << value << " ";}
    std::cout << "\nSim Runtime: " << simTimer*dt << " seconds";
    return 0;
}
//To-Do: Get motion to work. All data transfers are working, but model is not moving.