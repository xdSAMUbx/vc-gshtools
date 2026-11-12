// tide_system.cpp: parsing of permanent-tide system names and their canonical spelling.
#include "tide_system.hpp"
#include "enum_parse.hpp"

namespace gsh
{
  // Alias table: the canonical names are marked `true` and are the only ones listed in error
  // messages. Aliases cover the Spanish names used in older data files and other spelling variants.
  TideSystem parse_tide(std::string_view s)
  {
    return detail::parse_enum<TideSystem>(s, "Tide system",
                                          {{"tide_free", TideSystem::TideFree, true},
                                           {"sin_marea", TideSystem::TideFree, false},
                                           {"libre_de_marea", TideSystem::TideFree, false},
                                           {"tidefree", TideSystem::TideFree, false},
                                           {"non_tidal", TideSystem::TideFree, false},
                                           {"nontidal", TideSystem::TideFree, false},
                                           {"zero_tide", TideSystem::ZeroTide, true},
                                           {"marea_cero", TideSystem::ZeroTide, false},
                                           {"zero", TideSystem::ZeroTide, false},
                                           {"mean_tide", TideSystem::MeanTide, true},
                                           {"marea_media", TideSystem::MeanTide, false},
                                           {"mean", TideSystem::MeanTide, false}});
  }

  const char *to_string(TideSystem v)
  {
    switch (v)
    {
    case TideSystem::TideFree: return "tide_free";
    case TideSystem::ZeroTide: return "zero_tide";
    case TideSystem::MeanTide: return "mean_tide";
    }
    return "?"; // unreachable for valid enum values; kept for defensive builds
  }
}
