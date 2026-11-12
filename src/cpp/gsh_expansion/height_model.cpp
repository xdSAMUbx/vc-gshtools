// height_model.cpp: height source from HC/HS coefficients (moved here from gsh_utils/height.cpp).
// It uses the point synthesis of shs.hpp with q = 1 and no GM/r factor, so the output is the height itself.
#include "height_model.hpp"
#include <algorithm>
#include <cstdint>
#include <vector>
#include "shs.hpp"

namespace gsh
{
  namespace
  {
    // Height from HC/HS coefficients, with q = 1 and no GM/r factor.
    class ModelSource final : public HeightSource
    {
    public:
      ModelSource(Model m, int nmax)
      {
        // The topography model is not a gravity model: no GM scaling. Setting R = 1 and evaluating at r = 1
        // gives q = R / r = 1 in the kernels, so the raw sum S_0 is the height itself.
        m.R = 1.0;
        prep_ = prepare(m, 0, nmax < 0 ? m.nmax : nmax);
      }
      // The geodetic-latitude argument is ignored: the model is evaluated at the spherical latitude of pts.
      void evaluate(const Points &pts, const AlignedVec<double> &, std::vector<double> &H, const SynthOptions &o) const override
      {
        // Same angles, radius set to 1 so that q = 1 (see the constructor).
        Points unit = pts;
        std::fill(unit.r.begin(), unit.r.end(), 1.0);
        const auto raw = synth_points(prep_, unit, 0, o);
        H.assign(raw.f[0].begin(), raw.f[0].end());
      }

    private:
      Prepared prep_;
    };
  }

  std::shared_ptr<HeightSource> height_from_model(Model model, int nmax)
  {
    return std::make_shared<ModelSource>(std::move(model), nmax);
  }
}
