#pragma once
#include <string_view>
#include <array>

class Model
{
public:
  const std::string_view Name;

  const double Ixx; // moment of inertia, roll axis
  const double Iyy; // moment of inertia, pitch axis
  const double Izz; // moment of inertia, yaw axis
  const double Ixy; // product of inertia
  const double Iyz; // product of inertia
  const double Izx; // product of inertia

  const double m;    // mass
  const double xbar; // center of gravity offset, x
  const double ybar; // center of gravity offset, y
  const double zbar; // center of gravity offset, z
  const double S;      // reference area
  const double b;      // reference span
  const double cbar;   // reference chord length
  const double l_body; // body length
  const double l_nose; // nose length
  const double d;      // diameter
  const double CL;   // lift coefficient
  const double CD;   // drag coefficient
  const double CY;   // side force coefficient
  const double Cl;   // rolling moment coefficient
  const double Cm;   // pitching moment coefficient
  const double Cn;   // yawing moment coefficient
  const double Clp;  // roll damping derivative
  const double Clr;  // roll-due-to-yaw-rate derivative
  const double Cmq;  // pitch damping derivative
  const double Cnp;  // yaw-due-to-roll-rate derivative
  const double Cnr;  // yaw damping derivative

  const std::array<double, 3> length;
  ~Model() = default;
protected:
  Model(
    std::string_view Name,
    double Ixx, double Iyy, double Izz,
    double Ixy, double Iyz, double Izx,
    double m, double xbar, double ybar, double zbar,
    double S, double b, double cbar,
    double l_body, double l_nose, double d,
    double CL, double CD, double CY,
    double Cl, double Cm, double Cn,
    double Clp, double Clr, double Cmq, double Cnp, double Cnr,
    std::array<double, 3> length
  )
    : Name{Name}
    , Ixx{Ixx}, Iyy{Iyy}, Izz{Izz}
    , Ixy{Ixy}, Iyz{Iyz}, Izx{Izx}
    , m{m}, xbar{xbar}, ybar{ybar}, zbar{zbar}
    , S{S}, b{b}, cbar{cbar}
    , l_body{l_body}, l_nose{l_nose}, d{d}
    , CL{CL}, CD{CD}, CY{CY}
    , Cl{Cl}, Cm{Cm}, Cn{Cn}
    , Clp{Clp}, Clr{Clr}, Cmq{Cmq}, Cnp{Cnp}, Cnr{Cnr}
    , length{length}
  {}
};
