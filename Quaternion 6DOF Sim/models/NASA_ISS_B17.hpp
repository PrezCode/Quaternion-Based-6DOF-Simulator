//NASA ISS for orbital check cases; NESC-RP-12-00770, V1.0
#pragma once
#include <iostream>
#include <cmath>
#include <array>
#include <string>
enum Axis{x = 0, y = 1, z = 2};
class InternationalSpaceStation{
    public:
    InternationalSpaceStation(){
        modelData[0] = I[x][x] = 1.02E8;
        modelData[1] = I[y][y] = 0.91E8;
        modelData[2] = I[z][z] = 1.64E8;
        modelData[3] = I[x][y] = 6.96E6;
        modelData[4] = I[y][z] = -5.90E5;
        modelData[5] = I[z][x] = 5.48E6;        
        modelData[6] = mass = 4E5;
        modelData[7] = cgOffset[x] = -3.0;
        modelData[8] = cgOffset[y] = -1.5;
        modelData[9] = cgOffset[z] = 4.0;
        modelData[10] = Aref = 0.0;
        modelData[11] = span = 0.0;
        modelData[12] = chord = 0.0;
        modelData[13] = length_body = 0.0;
        modelData[14] = length_nose = 0.0;
        modelData[15] = diameter = 0.0;
        modelData[16] = C_lift = 0.0;
        modelData[17] = C_drag = 0.0;
        modelData[18] = C_sideforce = 0.0;
        modelData[19] = C_l = 0.0;
        modelData[20] = C_m = 0.0;
        modelData[21] = C_n = 0.0;
        modelData[22] = C_lp = 0.0;
        modelData[23] = C_lr = 0.0;
        modelData[24] = C_mq = 0.0;
        modelData[25] = C_np = 0.0;
        modelData[26] = C_nr = 0.0;
    }
    void getModel(std::array<double, 27>& dataTable){
        for(int i = 0; i < dataTable.size(); ++i){dataTable[i] = modelData[i];}
    }
    private:
    std::array<double, 24> modelData{0.0};
    double slug2kg{14.5939}, sqft2sqm{0.092903}, slugftsq2kgmsq{1.35581795}, m2ft{3.28084}, ft2m{0.3048},
    I[3][3]{0.0}, mass{0.0}, Aref{0.0}, span{0.0}, chord{0.0}, cgOffset[3]{0.0}, length_body{0.0}, length_nose{0.0}, diameter{0.0},
    C_lift{0.0}, C_drag{0.0}, C_sideforce{0.0}, C_l{0.0}, C_m{0.0}, C_n{0.0}, C_lp{0.0}, C_lr{0.0},
    C_mq{0.0}, C_np{0.0}, C_nr{0.0};
};