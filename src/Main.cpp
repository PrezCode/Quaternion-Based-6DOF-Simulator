//Written by Prez, https://github.com/PrezCode
//All code created without the use of AI
//Main file to test quaternion based simulator.
//Assume all values should be converted to metric and radians
#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include <iomanip>
#include <string>
#include <chrono>
#include <filesystem>
#include "computers/QuaternionSimulation.hpp"
#include "models/NASA_ISS_B17.hpp"
#include "models/NASA_Brick_B13.hpp"
#include <numbers>

void printInitialState(std::array<double, 12> state){
        std::cout << std::endl;
        std::cout << std::setw(10) << " " 
        << std::setw(10) << "Xe" << std::setw(10) << "Ye" << std::setw(10) << "Ze" << std::setw(10) << "Ue" << std::setw(10) << "Ve" << std::setw(10) << "We" 
        << std::setw(10) << "P" << std::setw(10) << "Q" << std::setw(10) << "R"<< std::setw(10) << "Phi" << std::setw(10) << "Theta" << std::setw(10) << "Psi" 
        << std::endl;
        std::cout << std::setw(10) << "Initial: ";
        for (auto& i : state)
        {
            std::cout << std::setw(10) << i << " ";
        }
        std::cout << std::endl;
}

int main(){
    auto start = std::chrono::system_clock::now();
    std::filesystem::create_directories("matlab");
    std::ofstream dataOut("matlab/Object_Simulation_Data.txt");
    std::array<double, 12> finalState{0.0};
    double dt {0.01}, r2d {180.0 / std::numbers::pi}, d2r {std::numbers::pi / 180.0}, m2ft {3.28084}, ft2m {0.3048};
    std::array<double, 12>  initialState = {-4315967.74, 960356.20, 5167269.53, 129.091037, -7491.513855, 1452.515654, 0*d2r, 0.0*d2r, 0*d2r, 0*d2r, 0.0*d2r, 0*d2r};//Xe, Ye, Ze, Ue, Ve, We, P, Q, R, Phi, Theta, Psi; if using latlong, XY = LatLong, Z = Alt above sea level
    int simEnd{5600}, simTimer{0};
    //Create and upload model
    QuaternionSimulator testObject(Brick(), ObjectType::Miscellaneous, Earth::Spherical, CoordinateInput::Cartesian);
    testObject.SetState(initialState);
    //Display initial state
    printInitialState(initialState);
    //Run simulation
    testObject.initialization();
    while(simTimer < simEnd/dt){
        testObject.iterate(dt);
        testObject.GetState(finalState);
        dataOut << simTimer*dt << " ";
        for(int i = 0; i < 12; ++i){dataOut << finalState[i] << " ";}
        dataOut << std::endl;
        if(abs(finalState[2]) < 0.0){break;}    //to stop sim if it hits surface on flat earth
        //if(sqrt(finalState[0]*finalState[0] + finalState[1]*finalState[1] + finalState[2]*finalState[2]) < 6371008.2){break;} //To stop sim if it hits spherical earth surface
        simTimer++;
    }
    //Print final state
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> runtime = end - start;
    std::cout << std::setw(10) << "Final: ";
    for(int value : finalState){std::cout << std::setw(10) << value << " ";}
    std::cout << "\nSim Runtime: " << simTimer*dt << " seconds" 
        << "\nReal Runtime: " << runtime.count() << std::endl;
    return 0;
}
