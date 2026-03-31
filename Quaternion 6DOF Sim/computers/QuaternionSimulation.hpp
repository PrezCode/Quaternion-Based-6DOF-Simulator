//Written by Prez, https://github.com/PrezCode
//Formulas from https://apps.dtic.mil/sti/tr/pdf/ADA417259.pdf and "Missile Guidance and Control Systems" by George M. Siouris
//All code created without the use of AI
#pragma once
#include <iostream>
#include <cmath>
#include <string>
#include <array>
#include "Atmosphere.hpp"
enum Variables{X = 0, Y = 1, Z = 2, U = 0, V = 1, W = 2, P = 0, Q = 1, R = 2, L = 0, M = 1, N = 2};
enum ObjectType{Missile, Aircraft, Miscellaneous};
enum Earth{Flat, Spherical};
enum CoordinateInput{Cartesian, LatLong};
class QuaternionSimulator{
    public:
    QuaternionSimulator(std::array<double, 27>& model, ObjectType object, Earth shape, CoordinateInput coordType){
        std::cout << "Quaternion 6DOF Simulation:" << std::endl;
        if(object == Missile || object == Aircraft){isMiscellanous = false;}
            else if(object == Miscellaneous){isMiscellanous = true;}
        if(shape == Flat){flatEarth = true; std::cout << "Flat Earth ON" << std::endl;}
            else if(shape == Spherical){flatEarth = false; std::cout << "Flat Earth OFF" << std::endl;}
        if(coordType == Cartesian){isCartesian = true; std::cout << "Cartesian ON" << std::endl;}
            else if(coordType == LatLong){isCartesian = false; std::cout << "Cartesian OFF" << std::endl;}
        I[X][X] = model[0];
        I[Y][Y] = model[1];
        I[Z][Z] = model[2];
        I[X][Y] = model[3];
        I[Y][Z] = model[4];
        I[Z][X] = model[5];
        mass = model[6];
        cgOffset[X] = model[7];
        cgOffset[Y] = model[8];
        cgOffset[Z] = model[9];
        A_ref = model[10];
        span = model[11];
        chord = model[12];
        length_body = model[13];
        length_nose = model[14];
        diameter = model[15];
        C_lift = model[16];
        C_drag = model[17];
        C_sideforce = model[18];
        C_l = model[19];
        C_m = model[20];
        C_n = model[21];
        C_lp = model[22];
        C_lr = model[23];
        C_mq = model[24];
        C_np = model[25];
        C_nr = model[26];
    }
    void setState(double* initials){
        for(int i = 0; i < 12; ++i){state[i] = initials[i];}
        if(flatEarth == true){state[2] = -state[2];}    //forces NED coordinates
            else if((flatEarth == false) && (isCartesian == true)){convertCartesionToECEF();}
            else if((flatEarth == false) && (isCartesian == false)){
                for(int i = 0; i < 3; ++i){latlong[i] = state[i];}
                convertECEFtoCartesion();
                for(int i = 0; i < 3; ++i){state[i] = cartesian[i];}
            }
    }
    void getState(std::array<double, 12>& outStates){
        for(int i = 0; i < 12; ++i){outStates[i] = state[i];}
        if(flatEarth == true){outStates[2] = -outStates[2];}
        outStates[3] = V_e[X];
        outStates[4] = V_e[Y];
        outStates[5] = V_e[Z];
    }
    void initialization(){
        Vmag_e = sqrt(V_e[U]*V_e[U] + V_e[V]*V_e[V] + V_e[W]*V_e[W]);
        V_e[U] = Vmag_e*cos(theta)*cos(psi);
        V_e[V] = Vmag_e*cos(phi)*sin(psi);
        V_e[W] = -Vmag_e*sin(theta);
        generateQuaternions();
        generateDCM();
        V_b[U] = C_eb[0][0]*V_e[U] + C_eb[1][0]*V_e[V] + C_eb[2][0]*V_e[W];
        V_b[V] = C_eb[0][1]*V_e[U] + C_eb[1][1]*V_e[V] + C_eb[2][1]*V_e[W];
        V_b[W] = C_eb[0][2]*V_e[U] + C_eb[1][2]*V_e[V] + C_eb[2][2]*V_e[W];
        Vmag_b = sqrt(V_b[U]*V_b[U] + V_b[V]*V_b[V] + V_b[W]*V_b[W]);
    }
    void iterate(double dt){
        RK4Integrator(dt);
        unpackIteration();
    }
    
    private:
    void RK4Integrator(double h){
        double k[4][12]{0.0}, y0[12]{0.0};
        for(int i = 0; i < 12; ++i){y0[i] = state[i];}
        unpackIteration();
        generateInstance(h);
        for(int i = 0; i < 12; ++i){
            k[0][i] = derivative[i];
            state[i] = y0[i] + 0.5*h*k[0][i];
        }
        unpackIteration();
        generateInstance(h);
        for(int i = 0; i < 12; ++i){
            k[1][i] = derivative[i];
            state[i] = y0[i] + 0.5*h*k[1][i];
        }
        unpackIteration();
        generateInstance(h);
        for(int i = 0; i < 12; ++i){
            k[2][i] = derivative[i];
            state[i] = y0[i] + h*k[2][i];
        }
        unpackIteration();
        generateInstance(h);
        for(int i = 0; i < 12; ++i){
            k[3][i] = derivative[i];
        }
        for(int i = 0; i < 12; ++i){state[i] = y0[i] + h*(k[0][i] + 2.0*k[1][i] + 2.0*k[2][i] + k[3][i])/6.0;}
        //std::cout << std::setw(10) << " ";
        //for(int i = 0; i < 12; ++i){std::cout << std::setw(7) << state[i];}
        //std::cout << std::endl;
    }
    void generateInstance(double dt){
        fluid.generateComplexAtmosphere(state[2]);
        generateForces();
        generateStateDerivatives();
        processing(dt);
    }
    void processing(double dt){
        generateEarthVelocities();
        generateAngularDeltas(dt);
        generateMaclaurinSeries();
        generateQuaternionDerivatives();
        generateNormalizedQuaternions();
        generateDCM();
        generateAirData();
    }
    void generateQuaternions(){
        quat[0] = sin(psi/2.0)*sin(theta/2.0)*cos(phi/2.0) - cos(psi/2.0)*cos(theta/2.0)*sin(phi/2.0);
        quat[1] = -1.0*cos(psi/2.0)*sin(theta/2.0)*cos(phi/2.0) - sin(psi/2.0)*cos(theta/2.0)*sin(phi/2.0);
        quat[2] = -1.0*sin(psi/2.0)*cos(theta/2.0)*cos(phi/2.0) + cos(psi/2.0)*sin(theta/2.0)*sin(phi/2.0);
        quat[3] = -1.0*cos(psi/2.0)*cos(theta/2.0)*cos(phi/2.0) - sin(psi/2.0)*sin(theta/2.0)*sin(phi/2.0);
        //check if quaternions are working
        //double quatCheck = quat[0]*quat[0] + quat[1]*quat[1] + quat[2]*quat[2] + quat[3]*quat[3];
        //std::cout << "Initialized quaternion magnitude sum: " << quatCheck << std::endl;
    }
    void generateForces(){
        if(flatEarth == false){g = -generateGravity();}
            else if(flatEarth == true){g = gStandard;}
        double qBar = 0.5*fluid.getAirDensity()*Vmag_b*Vmag_b;
        if(isMiscellanous == false){generateMachDrag(Vmag_b/fluid.getSoundBarrier(), qBar);}
        lift = C_lift*qBar*A_ref;
        drag = C_drag*qBar*A_ref;
        sideforce = C_sideforce*qBar*A_ref;
        //Force applied to the body axis
        force[X] = -(C_wb[0][0]*drag + C_wb[0][1]*sideforce + C_wb[0][2]*lift) + thrust;
        force[Y] = -(C_wb[1][0]*drag + C_wb[1][1]*sideforce + C_wb[1][2]*lift);
        force[Z] = -(C_wb[2][0]*drag + C_wb[2][1]*sideforce + C_wb[2][2]*lift);
        //Moments applied to the body axis
        if(Vmag_b > 0.0){
            C_l = ((C_lp*omega[P] + C_lr*omega[R])*span)/(2.0*Vmag_b);
            C_m = (C_mq*omega[Q]*chord)/(2.0*Vmag_b);
            C_n = ((C_np*omega[P] + C_nr*omega[R])*span)/(2.0*Vmag_b);
        }
        moment[L] = C_l*qBar*span*A_ref;
        moment[M] = C_m*qBar*chord*A_ref;
        moment[N] = C_n*qBar*span*A_ref;
    }
    void unpackIteration(){
        for(int i = 0; i < 3; ++i){
            cartesian[i] = state[i];
            V_b[i] = state[i+3];
            omega[i] = state[i+6];
        }
        phi = state[9];
        theta = state[10];
        psi = state[11];
        if(flatEarth == false){convertCartesionToECEF();}
    }
    void generateAirData(){
        Vmag_b = sqrt(V_b[U]*V_b[U] + V_b[V]*V_b[V] + V_b[W]*V_b[W]);
        Vmag_e = sqrt(V_e[U]*V_e[U] + V_e[V]*V_e[V] + V_e[W]*V_e[W]);
        if(Vmag_b == 0.0){
            alpha = 0.0;
            beta = 0.0;
        }else{
            alpha = atan2(V_b[W], V_b[U]);
            beta = asin(V_b[V]/Vmag_b);
        }
        gamma = atan2(-1.0*V_e[W], sqrt(V_e[U]*V_e[U] + V_e[V]*V_e[V]));
    }
    void generateStateDerivatives(){
        if(flatEarth == true){
            //Displacement Velocity
            derivative[0] = C_bn[0][0]*V_b[U] + C_bn[1][0]*V_b[V] + C_bn[2][0]*V_b[W];
            derivative[1] = C_bn[0][1]*V_b[U] + C_bn[1][1]*V_b[V] + C_bn[2][1]*V_b[W];
            derivative[2] = C_bn[0][2]*V_b[U] + C_bn[1][2]*V_b[V] + C_bn[2][2]*V_b[W];
            //Body Accelerations
            derivative[3] = force[X]/mass + C_eb[2][0]*g - omega[Q]*V_b[W] + omega[R]*V_b[V];
            derivative[4] = force[Y]/mass + C_eb[2][1]*g - omega[R]*V_b[U] + omega[P]*V_b[W];
            derivative[5] = force[Z]/mass + C_eb[2][2]*g - omega[P]*V_b[V] + omega[Q]*V_b[U];
        }else if(flatEarth == false){
            double mu = G_constant*(pow(massEarth, 3)/pow(mass + massEarth, 2));
            //Displacement Velocity
            derivative[0] = C_eb[0][0]*V_b[U] + C_eb[1][0]*V_b[V] + C_eb[2][0]*V_b[W];
            derivative[1] = C_eb[0][1]*V_b[U] + C_eb[1][1]*V_b[V] + C_eb[2][1]*V_b[W];
            derivative[2] = C_eb[0][2]*V_b[U] + C_eb[1][2]*V_b[V] + C_eb[2][2]*V_b[W];
            //Body Accelerations
            derivative[3] = -mu/pow(cartesian[X]*cartesian[X] + cartesian[Y]*cartesian[Y] + cartesian[Z]*cartesian[Z], 1.5)*cartesian[X];
            derivative[4] = -mu/pow(cartesian[X]*cartesian[X] + cartesian[Y]*cartesian[Y] + cartesian[Z]*cartesian[Z], 1.5)*cartesian[Y];
            derivative[5] = -mu/pow(cartesian[X]*cartesian[X] + cartesian[Y]*cartesian[Y] + cartesian[Z]*cartesian[Z], 1.5)*cartesian[Z];
        }   
        //Rotation Accelerations
        derivative[6] = (I[X][Z]*(I[X][X] - I[Y][Y] + I[Z][Z])*omega[P]*omega[Q] - (I[Z][Z]*(I[Z][Z] - I[Y][Y]) + I[X][Z]*I[X][Z])*omega[Q]*omega[R] + I[Z][Z]*moment[L] + I[X][Z]*moment[N])/(I[X][X]*I[Z][Z] - I[X][Z]*I[X][Z]);
        derivative[7] = ((I[Z][Z] - I[X][X])*omega[R]*omega[P] - I[X][Z]*(omega[P]*omega[P] - omega[R]*omega[R]) + moment[M])/I[Y][Y];
        derivative[8] = (-I[X][Z]*(I[X][X] - I[Y][Y] + I[Z][Z])*omega[Q]*omega[R] + (I[X][X]*(I[X][X] - I[Y][Y]) + I[X][Z]*I[X][Z])*omega[P]*omega[Q] + I[X][Z]*moment[L] + I[X][X]*moment[N])/(I[X][X]*I[Z][Z] - I[X][Z]*I[X][Z]);
        //Euler Angle Velocities
        derivative[9] = (omega[P] + sin(phi)*tan(theta))*omega[Q] + cos(phi)*tan(theta)*omega[R];
        derivative[10] = cos(phi)*omega[Q] - sin(phi)*omega[R];
        derivative[11] = (sin(phi)/cos(theta))*omega[Q] + (cos(phi)/cos(theta))*omega[R];
    }
    void generateEarthVelocities(){
        V_e[X] = C_bn[0][0]*V_b[U] + C_bn[1][0]*V_b[V] + C_bn[2][0]*V_b[W];
        V_e[Y] = C_bn[0][1]*V_b[U] + C_bn[1][1]*V_b[V] + C_bn[2][1]*V_b[W];
        V_e[Z] = C_bn[0][2]*V_b[U] + C_bn[1][2]*V_b[V] + C_bn[2][2]*V_b[W];
    }
    void generateAngularDeltas(double dt){
        deltaPhi  = omega[P]*dt;
        deltaTheta = omega[Q]*dt;
        deltaPsi = omega[R]*dt;
    }
    void generateQuaternionDerivatives(){
        generateMaclaurinSeries();
        quatDeriv[0] = quat[0]*Cn + quat[1]*Sn*deltaPsi + -1.0*quat[2]*Sn*deltaTheta + quat[3]*Sn*deltaPhi;
        quatDeriv[1] = -1.0*quat[0]*Sn*deltaPsi + quat[1]*Cn + quat[2]*Sn*deltaPhi + quat[3]*Sn*deltaTheta;
        quatDeriv[2] = quat[0]*Sn*deltaTheta + -1.0*quat[1]*Sn*deltaPhi + quat[2]*Cn + quat[3]*Sn*deltaPsi;
        quatDeriv[3] = -1.0*quat[0]*Sn*deltaPhi + -1.0*quat[1]*Sn*deltaTheta + -1.0*quat[2]*Sn*deltaPsi;
    }
    void generateMaclaurinSeries(){
        Cn = 1.0 - (deltaPhi*deltaPhi + deltaTheta*deltaTheta + deltaPsi*deltaPsi)/8.0 
            + (powl(deltaPhi, 4) + powl(deltaTheta, 4) + powl(deltaPsi, 4))/384;
        Sn = 0.5 - (deltaPhi*deltaPhi + deltaTheta*deltaTheta + deltaPsi*deltaPsi)/48.0;
    }
    void generateNormalizedQuaternions(){
        double normalizer = 0.5*(3.0 - quat[0]*quat[0] - quat[1]*quat[1] - quat[2]*quat[2] - quat[3]*quat[3]);
        for(int i = 0; i < 4; ++i){quat[i] *= normalizer;}
    } 
    void generateDCM(){
        //earth-to-body
        C_eb[0][0] = quat[0]*quat[0] - quat[1]*quat[1] - quat[2]*quat[2] + quat[3]*quat[3];
        C_eb[0][1] = 2*(quat[0]*quat[1] - quat[2]*quat[3]);
        C_eb[0][2] = 2*(quat[0]*quat[2] + quat[1]*quat[3]);
        C_eb[1][0] = 2*(quat[0]*quat[1] + quat[2]*quat[3]);
        C_eb[1][1] = -1.0*quat[0]*quat[0] + quat[1]*quat[1] - quat[2]*quat[2] + quat[3]*quat[3];
        C_eb[1][2] = 2*(quat[1]*quat[2] - quat[0]*quat[3]);
        C_eb[2][0] = 2*(quat[0]*quat[2] - quat[1]*quat[3]);
        C_eb[2][1] = 2*(quat[1]*quat[2] + quat[0]*quat[3]);
        C_eb[2][2] = -1.0*quat[0]*quat[0] - quat[1]*quat[1] + quat[2]*quat[2] + quat[3]*quat[3];
        //body-to-navigation
        C_bn[0][0] = quat[0]*quat[0] + quat[1]*quat[1] - quat[2]*quat[2] - quat[3]*quat[3];
        C_bn[0][1] = 2*(quat[0]*quat[1] - quat[2]*quat[3]);
        C_bn[0][2] = 2*(quat[0]*quat[2] + quat[1]*quat[3]);
        C_bn[1][0] = 2*(quat[0]*quat[1] + quat[2]*quat[3]);
        C_bn[1][1] = quat[0]*quat[0] - quat[1]*quat[1] + quat[2]*quat[2] - quat[3]*quat[3];
        C_bn[1][2] = 2*(quat[1]*quat[2] - quat[0]*quat[3]);
        C_bn[2][0] = 2*(quat[0]*quat[2] - quat[1]*quat[3]);
        C_bn[2][1] = 2*(quat[1]*quat[2] + quat[0]*quat[3]);
        C_bn[2][2] = quat[0]*quat[0] - quat[1]*quat[1] - quat[2]*quat[2] + quat[3]*quat[3];
        //Wind-to-Body
        C_wb[0][0] = cos(alpha)*cos(beta);
        C_wb[0][1] = -cos(alpha)*sin(beta);
        C_wb[0][2] = -sin(alpha);
        C_wb[1][0] = sin(beta);
        C_wb[1][1] = cos(beta);
        C_wb[1][2] = 0;
        C_wb[2][0] = sin(alpha)*cos(beta);
        C_wb[2][1] = -sin(alpha)*sin(beta);
        C_wb[2][2] = cos(alpha);
    }
    double wrapToPi(double angle){//Force angle to be -pi < x < pi
        while(angle > M_PI){angle = angle - M_PI;}
        while(angle < M_PI){angle = angle + M_PI;}
        return angle;
    }
    double wrapTo2Pi(double angle){//Force angle to be -2pi < x < 2pi
        while(angle > 2*M_PI){angle = angle - 2*M_PI;}
        while(angle < 2*M_PI){angle = angle + 2*M_PI;}
        return angle;
    }
    void generateMachDrag(double mach, double qBar){
        double C_dragWave, C_dragCoast, C_dragPowered, C_dragFriction;
        if(mach > 1){
            C_dragWave = 3.6/((length_nose/diameter)*(mach - 1.0) + 3.0);
            if(thrust > 0.0){C_dragCoast = 0.25/mach;}
                else{C_dragPowered = (1.0 - A_nozzleExit/A_ref)*(0.25/mach);}
        }else if(mach <= 1){
            C_dragWave = 0.0;
            if(thrust > 0.0){C_dragPowered = (1.0 - A_nozzleExit/A_ref)*(0.12 + 0.13*mach*mach);}
                else{C_dragCoast = (0.12 + 0.13*mach*mach);}
        }
        C_dragFriction = 0.053*(length_body/diameter)*pow(mach/(qBar*length_body), 0.2);
        C_drag = C_dragWave + C_dragCoast + C_dragPowered + C_dragFriction;
    }
    void convertECEFtoCartesion(){
        cartesian[X] = (latlong[Z] + radiusEarth)*cos(latlong[X])*cos(latlong[Y]);
        cartesian[Y] = (latlong[Z] + radiusEarth)*cos(latlong[X])*sin(latlong[Y]);
        cartesian[Z] = (latlong[Z] + radiusEarth)*sin(latlong[X]);
    }
    void convertCartesionToECEF(){
        latlong[X] = asin(cartesian[Z]/sqrtl(cartesian[X]*cartesian[X] + cartesian[Y]*cartesian[Y] + cartesian[Z]*cartesian[Z]));
        latlong[Y] = atan2(cartesian[Y], cartesian[X]);
        latlong[Z] = sqrtl(cartesian[X]*cartesian[X] + cartesian[Y]*cartesian[Y] + cartesian[Z]*cartesian[Z]) - radiusEarth;
    }
    double generateGravity(){return (G_constant*mass*massEarth)/pow(latlong[Z] + radiusEarth, 2.0);}
    Atmosphere fluid;
    std::array<double, 12> state{0.0};
    bool isMiscellanous, flatEarth, isCartesian;
    double quat[4]{0.0}, quatDeriv[4]{0.0}, mass{0.0}, massEarth{5.9722E24/*kg*/}, radiusEarth{6371008.2/*m*/}, earthRotation{7292115.0E-11/*rad/s*/}, derivative[12]{0.0}, 
    phi{0.0}, theta{0.0}, psi{0.0}, deltaPhi{0.0}, deltaTheta{0.0}, deltaPsi{0.0}, Cn{0.0}, Sn{0.0}, cartesian[3]{0.0}, latlong[3]{0.0},
    C_eb[3][3]{0.0}, C_bn[3][3]{0.0}, C_wb[3][3]{0.0}, C_lift{0.0}, C_drag{0.0}, C_sideforce{0.0},  
    V_b[3]{0.0}, Vmag_b{0.0}, V_e[3]{0.0}, Vmag_e{0.0},
    omega[3]{0.0}, beta{0.0}, alpha{0.0}, gamma{0.0},
    force[3]{0.0}, lift{0.0}, drag{0.0}, sideforce{0.0}, thrust{0.0}, moment[3]{0.0},
    I[3][3]{0.0}, gStandard{9.80665/*m/s^2*/}, g{0.0}, G_constant{6.674E-11/*m^3/kgs^2*/}, C_l{0.0}, C_m{0.0}, C_n{0.0}, C_lp{0.0}, C_lr{0.0}, C_mq{0.0}, C_np{0.0}, C_nr{0.0},
    A_ref{0.0}, span{0.0}, chord{0.0}, cgOffset[3]{0.0}, length_body{0.0}, length_nose{0.0}, diameter{0.0}, A_nozzleExit{0.0};
};
