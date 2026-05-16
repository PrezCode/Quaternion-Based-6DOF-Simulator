#pragma once

namespace Units
{
  enum Axis 
  {
    x, y, z
  };

  inline constexpr double ft2m = 0.3048;
  inline constexpr double m2ft = 1.0 / ft2m;
  inline constexpr double sqft2sqm = ft2m * ft2m;
  inline constexpr double in2m = 0.0254;
  inline constexpr double m2in = 1.0 / in2m;
  inline constexpr double slug2kg = 14.5939;
  inline constexpr double slugftsq2kgmsq = slug2kg * ft2m * ft2m;

  constexpr double FeetToMeters(double ft) { return ft * ft2m; }
  constexpr double MetersToFeet(double m) { return m * m2ft; }
  constexpr double InchesToMeters(double in) { return in * in2m; }
  constexpr double MetersToInches(double m) { return m * m2in; }
  constexpr double SqFeetToSqMeters(double sqft) { return sqft * sqft2sqm; }
  constexpr double SlugsToKg(double slugs) { return slugs * slug2kg; }
  constexpr double SlugFtSqToKgMSq(double slugftsq) { return slugftsq * slugftsq2kgmsq; }
}
