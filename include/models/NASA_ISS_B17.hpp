//NASA ISS for orbital check cases; NESC-RP-12-00770, V1.0
#pragma once
#include <array>

enum Axis{x = 0, y = 1, z = 2};
class InternationalSpaceStation{
    public:
    InternationalSpaceStation(){
        m_modelData[0] = m_I[x][x] = 1.02E8;
        m_modelData[1] = m_I[y][y] = 0.91E8;
        m_modelData[2] = m_I[z][z] = 1.64E8;
        m_modelData[3] = m_I[x][y] = 6.96E6;
        m_modelData[4] = m_I[y][z] = -5.90E5;
        m_modelData[5] = m_I[z][x] = 5.48E6;
        m_modelData[6] = m_mass = 4E5;
        m_modelData[7] = m_cgOffset[x] = -3.0;
        m_modelData[8] = m_cgOffset[y] = -1.5;
        m_modelData[9] = m_cgOffset[z] = 4.0;
        m_modelData[10] = m_Aref = 0.0;
        m_modelData[11] = m_span = 0.0;
        m_modelData[12] = m_chord = 0.0;
        m_modelData[13] = m_length_body = 0.0;
        m_modelData[14] = m_length_nose = 0.0;
        m_modelData[15] = m_diameter = 0.0;
        m_modelData[16] = m_C_lift = 0.0;
        m_modelData[17] = m_C_drag = 0.0;
        m_modelData[18] = m_C_sideforce = 0.0;
        m_modelData[19] = m_C_l = 0.0;
        m_modelData[20] = m_C_m = 0.0;
        m_modelData[21] = m_C_n = 0.0;
        m_modelData[22] = m_C_lp = 0.0;
        m_modelData[23] = m_C_lr = 0.0;
        m_modelData[24] = m_C_mq = 0.0;
        m_modelData[25] = m_C_np = 0.0;
        m_modelData[26] = m_C_nr = 0.0;
    }
    void GetModel(std::array<double, 27>& dataTable){
        for(int i = 0; i < dataTable.size(); ++i){dataTable[i] = m_modelData[i];}
    }
    private:
    std::array<double, 24> m_modelData{0.0};
    double m_I[3][3]{0.0};
    double m_mass{0.0};
    double m_Aref{0.0};
    double m_span{0.0};
    double m_chord{0.0};
    double m_cgOffset[3]{0.0};
    double m_length_body{0.0};
    double m_length_nose{0.0};
    double m_diameter{0.0};
    double m_C_lift{0.0};
    double m_C_drag{0.0};
    double m_C_sideforce{0.0};
    double m_C_l{0.0};
    double m_C_m{0.0};
    double m_C_n{0.0};
    double m_C_lp{0.0};
    double m_C_lr{0.0};
    double m_C_mq{0.0};
    double m_C_np{0.0};
    double m_C_nr{0.0};
};
