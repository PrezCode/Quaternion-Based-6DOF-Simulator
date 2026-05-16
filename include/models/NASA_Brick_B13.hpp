// NASA Brick; NESC-RP-12-00770, V1.0
#pragma once
#include "Model.hpp"
#include "units/Units.hpp"
#include <array>

class Brick : public Model
{
public:
  Brick() : 
  Model (
    "Brick",
    Units::SlugFtSqToKgMSq(0.001894220),
    Units::SlugFtSqToKgMSq(0.006211019),
    Units::SlugFtSqToKgMSq(0.007194665),
    Units::SlugFtSqToKgMSq(0.0),
    Units::SlugFtSqToKgMSq(0.0),
    Units::SlugFtSqToKgMSq(0.0),
    Units::SlugsToKg(0.155404754),
    Units::MetersToFeet(0.0),
    Units::MetersToFeet(0.0),
    Units::MetersToFeet(0.0),
    Units::SqFeetToSqMeters(0.22222),
    Units::FeetToMeters(0.33333),
    Units::FeetToMeters(0.66667),
    0.0, 0.0, 0.0,
    0.0, 0.01, 0.0, 0.0, 0.0, 0.0, -1.0,
    0.0, -1.0, 0.0, -1.0, 
    {Units::InchesToMeters(8), 
    Units::InchesToMeters(4), 
    Units::InchesToMeters(2.25)})
  {

  }
  ~Brick() = default;
};
