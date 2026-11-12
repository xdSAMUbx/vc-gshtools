// height_model.hpp: height source H(φ, λ) from a spherical-harmonic model with HC/HS coefficients (topography).
// It lives in gsh_expansion because it evaluates the synthesis of shs.hpp. The HeightSource interface and the
// regional grid source are in gsh_utils/height.hpp.
#pragma once
#include <memory>
#include "height.hpp"
#include "model.hpp"

namespace gsh
{
  // Wraps an HC/HS model as a HeightSource. The model is evaluated with r = R, so q = 1, and the output is H in metres
  // (the sum of HC cos(mλ) + HS sin(mλ) times P̄ of Bucha and Janák (2019)). nmax < 0 uses the whole model.
  // The source evaluates on the sphere: it uses the spherical latitude of each point, not the geodetic one.
  std::shared_ptr<HeightSource> height_from_model(Model model, int nmax = -1);
}
