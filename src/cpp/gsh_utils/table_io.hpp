// table_io.hpp: fast reader for delimited numeric tables (for example observation or gravity files).
//
// Behaviour in short: the file is read in one piece; the first `header_lines` lines are skipped; for
// each remaining line the requested columns are parsed as doubles. A line in which any requested
// field is not a number is recorded in Table::skipped_lines and is NOT repaired, interpolated or
// guessed. Blank lines are ignored and not reported.
#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "core.hpp"

namespace gsh
{
  // Columns of numeric values extracted from a table.
  struct Table
  {
    std::vector<std::vector<double>> columns;  // columns[k] holds the values of the k-th requested column
    std::vector<std::int64_t> skipped_lines;   // 1-based line numbers of rejected lines, in file order
  };

  // Reads the table at `path`.
  //   sep: field delimiter; '\0' means fields are separated by runs of spaces or tabs.
  //   header_lines: number of leading lines to skip (counted in lines, not records).
  //   cols: zero-based indices of the fields to keep, in output order. Fields are split up to the
  //         largest index; a line shorter than that index is rejected.
  // Numbers are parsed with std::from_chars; a single leading '+' is accepted. Trailing CR, space and
  // tab characters are trimmed from each line. Throws std::runtime_error if the file cannot be
  // opened or read completely.
  Table read_table(const std::string &path, char sep, int header_lines, const std::vector<int> &cols);
}
