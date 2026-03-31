#pragma once
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

class MatrixComputation{
    public:
    MatrixComputation(){}
    void multiply(std::vector<std::vector<double>>& A, std::vector<double>& B){
        size_t rows = sizeof(A)/sizeof(A[0]);
        size_t columns = sizeof(B)/sizeof(B[0]);
        std::cout << "Rows: " << rows << std::endl;
        std::cout << "columns: " << columns << std::endl;
        if(rows != columns){std::cout << "Rows and Columns are not the same" << std::endl;}
            else if(rows == columns){std::cout << "Rows and Columns are the same" << std::endl;}
    }
    private:
};