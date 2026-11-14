// table_io.cpp: implementation of the delimited numeric table reader declared in table_io.hpp.
#include "table_io.hpp"
#include <algorithm>
#include <charconv>
#include <cstdio>
#include <cstring>
#include <stdexcept>

namespace gsh
{
  // Parsing outline:
  //   - The whole file is loaded into memory (a single fread).
  //   - The buffer is scanned line by line with memchr. Line numbers are 1-based and include the
  //     header lines, so Table::skipped_lines reports the line number as seen in an editor.
  //   - For each line, fields are walked from left to right until the largest requested index.
  //     With sep == '\0' the separators are runs of spaces/tabs; otherwise every occurrence of sep
  //     separates two fields (so empty fields are possible and make the line invalid if requested).
  //   - A requested field is converted with std::from_chars. A single leading '+' is skipped, since
  //     from_chars does not accept it. If any requested field fails, the whole line is rejected.
  //   - Accepted lines append one value per requested column, so all columns stay the same length.
  Table read_table(const std::string &path, char sep, int header_lines, const std::vector<int> &cols)
  {
    std::FILE *f = std::fopen(path.c_str(), "rb");
    if (!f)
      throw std::runtime_error("Cannot open: " + path);
    const long long size = gsh::detail::stream_size(f);
    std::string buf(static_cast<std::size_t>(size), '\0');
    const std::size_t got = std::fread(buf.data(), 1, buf.size(), f);
    std::fclose(f);
    if (got != buf.size())
      throw std::runtime_error("Incomplete read: " + path);

    Table tab;
    tab.columns.resize(cols.size());
    int maxcol = 0;
    for (int c : cols)
      maxcol = std::max(maxcol, c);
    std::vector<double> row(cols.size()); // values of the line being parsed

    const char *p = buf.data(), *const end = p + buf.size();
    std::int64_t lineno = 0;
    while (p < end)
    {
      const char *eol = static_cast<const char *>(std::memchr(p, '\n', end - p));
      if (!eol)
        eol = end;
      ++lineno;
      const char *q = p;   // start of the current line
      p = eol + 1;
      if (lineno <= header_lines)
        continue;
      // Trim trailing CR, space and tab (Windows line endings and padding).
      const char *lend = eol;
      while (lend > q && (lend[-1] == '\r' || lend[-1] == ' ' || lend[-1] == '\t'))
        --lend;
      if (lend == q)
        continue; // blank line: neither a record nor an error

      bool ok = true;
      int field = 0;
      while (ok && field <= maxcol)
      {
        // Whitespace-delimited mode: skip leading blanks before the field.
        if (sep == '\0')
          while (q < lend && (*q == ' ' || *q == '\t'))
            ++q;
        const char *tok = q;
        if (sep == '\0')
          while (q < lend && *q != ' ' && *q != '\t')
            ++q;
        else
          while (q < lend && *q != sep)
            ++q;
        // An empty field at the end of the line means the line has fewer fields than requested.
        if (tok == q && q >= lend)
        {
          ok = false;
          break;
        }
        for (std::size_t c = 0; c < cols.size(); ++c)
          if (cols[c] == field)
          {
            const char *first = (*tok == '+') ? tok + 1 : tok;
            auto r = std::from_chars(first, q, row[c]);
            // The conversion must consume the whole field; otherwise the line is rejected.
            if (r.ec != std::errc() || r.ptr != q)
              ok = false;
          }
        // Move past the delimiter (only in explicit-separator mode, and only if one follows).
        if (sep != '\0' && q < lend)
          ++q;
        ++field;
      }
      if (!ok)
      {
        tab.skipped_lines.push_back(lineno);
        continue;
      }
      for (std::size_t c = 0; c < cols.size(); ++c)
        tab.columns[c].push_back(row[c]);
    }
    return tab;
  }
}
