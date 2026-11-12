// enum_parse.hpp: generic lookup table for parsing enumeration names from text. Matching is
// case-insensitive, treats '-' and ' ' as '_', and trims surrounding whitespace.
// Used by tide_system.cpp.
#pragma once
#include <cctype>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <string_view>

namespace gsh::detail
{
  // Normalizes an enumeration name: trims ' ', '\t', '\r', '\n' at both ends, converts ASCII
  // letters to lower case, and maps '-' and ' ' to '_'. Other characters are kept as is.
  // Returns an empty string for an all-whitespace input. Does not throw.
  inline std::string normalize(std::string_view s)
  {
    const std::size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string_view::npos)
      return {};
    const std::size_t e = s.find_last_not_of(" \t\r\n");
    std::string out;
    out.reserve(e - b + 1);
    for (char c : s.substr(b, e - b + 1))
    {
      if (c == '-' || c == ' ')
        out.push_back('_');
      else
        out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return out;
  }

  // One row of a parsing table. `canonical` marks the name that is listed in error messages;
  // other rows are aliases (for example Spanish names or other spellings).
  template <class E>
  struct Entry
  {
    std::string_view name; // already normalized (lower case, '_' separators)
    E value;
    bool canonical;
  };

  // Looks up `text` in `table` (after normalize) and returns the matching value.
  // `what` is the human-readable name of the enumeration, used in the error message.
  // Throws std::invalid_argument listing the canonical names if no row matches.
  template <class E>
  E parse_enum(std::string_view text, std::string_view what, std::initializer_list<Entry<E>> table)
  {
    const std::string key = normalize(text);
    for (const auto &e : table)
      if (key == e.name)
        return e.value;

    std::string valid;
    for (const auto &e : table)
    {
      if (!e.canonical)
        continue;
      if (!valid.empty())
        valid += ", ";
      valid += e.name;
    }
    throw std::invalid_argument(std::string(what) + " unknown: '" + std::string(text) +
                                "'. Valid: " + valid);
  }
} // namespace gsh::detail
