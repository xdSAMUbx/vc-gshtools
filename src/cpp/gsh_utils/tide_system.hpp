// tide_system.hpp: permanent-tide systems (tide-free, zero-tide, mean-tide) and their names.
// Lives in the gsh_utils base layer so that the tide conversions (tides.hpp) can use the same type.
//
// Background (Mäkinen 2021, section 2):
//   - TideFree  (non-tidal, NT): removes the direct luni-solar attraction and the permanent elastic
//               deformation. Standard for ITRF/WGS84 and GNSS coordinates.
//   - ZeroTide  (ZT): removes the direct attraction but keeps the permanent deformation.
//   - MeanTide  (MT): keeps both. Natural reference for mean sea level (altimetry, tide gauges).
#pragma once
#include <string_view>

namespace gsh
{
  // The three permanent-tide systems. Closed enumeration: switch statements over it need no default.
  enum class TideSystem { TideFree, ZeroTide, MeanTide };

  // Parses a system name. Accepted: the canonical names "tide_free", "zero_tide", "mean_tide" and
  // several aliases (for example "non_tidal", "zero", "mean"). Matching is case-insensitive and
  // accepts '-' or ' ' in place of '_'. Throws std::invalid_argument, listing the canonical names,
  // if the text matches no system.
  TideSystem parse_tide(std::string_view s);
  // Canonical name of a system: "tide_free", "zero_tide" or "mean_tide". Never returns nullptr.
  const char *to_string(TideSystem v);
}
