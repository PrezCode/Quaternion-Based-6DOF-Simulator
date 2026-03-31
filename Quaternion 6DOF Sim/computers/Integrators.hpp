#include <iostream>
#include <cmath>
#include <string>

class Integrators{
    public:
    Integrators(){}
    double SimpleEuler(double h, double y0, double dydt){
        double y = y0 + h*dydt;
        return y;
    }
    double ImprovedEuler(double h, double y0, double dydt_0, double dydt_1){
        double y = y0 + 0.5*h*(dydt_0 + dydt_1);
        return y;
    }
    double RK4(double h, double x, double y, double (*f)(double, double)){
        double k[4]{0.0}, y0{y};
        k[0] = f(x, y);
        k[1] = f(x + 0.5*h, y + 0.5*h*k[0]);
        k[2] = f(x + 0.5*h, y + 0.5*h*k[1]);
        k[3] = f(x + h, y + h*k[2]);
        double yn_1 = y0 + h*(k[0] + 2.0*k[1] + 2.0*k[2] + k[3])/6.0;
        return yn_1;
    }
    private:
};