// model.cpp: reader for ICGEM .gfc coefficient files and helpers for the packed coefficient model.
// See model.hpp for the format description and the error contract.
#include "model.hpp"
#include <algorithm>
#include <charconv>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string_view>

namespace gsh
{
  namespace
  {
    // Reads a whole file into memory with a single fread. Throws std::runtime_error if the file
    // cannot be opened or the read is short.
    std::string read_file(const std::string &path)
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
      return buf;
    }

    // Blank characters inside a data line. A '\n' is never blank here because lines are split first.
    inline bool is_blank(char c) { return c == ' ' || c == '\t' || c == '\r'; }

    // Parses a real number with optional leading '+' and optional exponent marker e/E or d/D (Fortran).
    // Skips leading blanks, advances p to the end of the token, and returns false if the token
    // is empty or not a number. The token is not modified on failure beyond the position update.
    inline bool parse_real(const char *&p, const char *end, double &out)
    {
      while (p < end && is_blank(*p))
        ++p;
      const char *tok = p;
      while (p < end && !is_blank(*p) && *p != '\n')
        ++p;
      if (tok == p)
        return false;
      const std::size_t len = static_cast<std::size_t>(p - tok);
      const char *first = (*tok == '+') ? tok + 1 : tok;
      auto r = std::from_chars(first, p, out);
      if (r.ec == std::errc() && r.ptr == p)
        return true;
      // Fallback for Fortran 'd'/'D' exponents: copy the token, replace the marker with 'e', retry.
      if (len >= 64)
        return false;
      char tmp[64];
      std::memcpy(tmp, tok, len);
      tmp[len] = '\0';
      for (char *q = tmp; *q; ++q)
        if (*q == 'd' || *q == 'D')
          *q = 'e';
      const char *f2 = (tmp[0] == '+') ? tmp + 1 : tmp;
      auto r2 = std::from_chars(f2, tmp + len, out);
      return r2.ec == std::errc() && r2.ptr == tmp + len;
    }

    // Parses an integer after skipping leading blanks; advances p past the digits.
    inline bool parse_int(const char *&p, const char *end, int &out)
    {
      while (p < end && is_blank(*p))
        ++p;
      auto r = std::from_chars(p, end, out);
      if (r.ec != std::errc())
        return false;
      p = r.ptr;
      return true;
    }

    // Second whitespace-separated token of a header line ("key value" -> "value"). Returns an empty
    // string if there is no such token.
    std::string header_value(std::string_view line)
    {
      std::size_t a = line.find_first_of(" \t");
      if (a == std::string_view::npos)
        return {};
      a = line.find_first_not_of(" \t", a);
      if (a == std::string_view::npos)
        return {};
      std::size_t b = line.find_first_of(" \t\r", a);
      return std::string(line.substr(a, b == std::string_view::npos ? b : b - a));
    }
  }

  // All-zero model with the given GM, R and degree. Both C and S are allocated with tri_size(nmax).
  Model Model::zeros(double GM, double R, int nmax)
  {
    Model m;
    m.GM = GM;
    m.R = R;
    m.nmax = m.file_nmax = nmax;
    m.C.assign(tri_size(nmax), 0.0);
    m.S.assign(tri_size(nmax), 0.0);
    return m;
  }

  // Copies the coefficients of degree <= new_nmax column by column. Metadata is copied unchanged,
  // except that nmax is set to new_nmax; file_nmax keeps the value declared in the source file.
  Model Model::truncated(int new_nmax) const
  {
    if (new_nmax > nmax)
      throw std::invalid_argument("truncated: degree larger than the model degree");
    Model out;
    out.GM = GM;
    out.R = R;
    out.file_nmax = file_nmax;
    out.name = name;
    out.norm = norm;
    out.tide_system = tide_system;
    out.nmax = new_nmax;
    out.C.resize(tri_size(new_nmax));
    out.S.resize(tri_size(new_nmax));
    for (int m = 0; m <= new_nmax; ++m)
    {
      // Column m holds n = m..nmax contiguously, so each column is a single block copy.
      const auto src = col_offset(m, nmax), dst = col_offset(m, new_nmax);
      const auto len = new_nmax - m + 1;
      std::copy_n(C.data() + src, len, out.C.data() + dst);
      std::copy_n(S.data() + src, len, out.S.data() + dst);
    }
    return out;
  }

  // Two-phase reader:
  //   1. Sequential header scan up to the "end_of_head" line; collects GM, R, max_degree and the
  //      optional text fields. Missing GM, R or max_degree (or a missing end_of_head) is an error.
  //   2. Parallel body scan. The body is split into chunks of about 4 MiB, each cut at a line start
  //      (a chunk that does not begin right after '\n' skips forward to the next line). Each thread
  //      writes only to the coefficient slots of its own lines, so no synchronization is needed.
  //      Lines that start with "gfc" and carry a malformed field are counted in `bad`; any such line
  //      makes the whole load fail after the parallel loop.
  namespace
  {
    Model load_gfc_impl(const std::string &path, int nmax, bool sigma);
  }
  Model load_gfc(const std::string &path, int nmax) { return load_gfc_impl(path, nmax, false); }
  Model load_gfc_sigma(const std::string &path, int nmax) { return load_gfc_impl(path, nmax, true); }

  namespace
  {
  Model load_gfc_impl(const std::string &path, int nmax, bool sigma)
  {
    const std::string buf = read_file(path);
    Model model;

    // Header, read sequentially up to end_of_head.
    std::size_t pos = 0;
    bool found_end = false;
    while (pos < buf.size())
    {
      std::size_t eol = buf.find('\n', pos);
      if (eol == std::string::npos)
        eol = buf.size();
      const std::string_view line(buf.data() + pos, eol - pos);
      pos = eol + 1;
      if (line.rfind("end_of_head", 0) == 0)
      {
        found_end = true;
        break;
      }
      // Parses the value of a numeric header key; throws on a malformed value.
      const auto parse_field = [&](double &dst)
      {
        const std::string v = header_value(line);
        const char *p = v.data();
        if (!parse_real(p, v.data() + v.size(), dst))
          throw std::runtime_error("Invalid GFC header line: " + std::string(line));
      };
      if (line.rfind("earth_gravity_constant", 0) == 0)
        parse_field(model.GM);
      else if (line.rfind("radius", 0) == 0)
        parse_field(model.R);
      else if (line.rfind("max_degree", 0) == 0)
      {
        double d = 0;
        parse_field(d);
        model.file_nmax = static_cast<int>(d);
      }
      else if (line.rfind("modelname", 0) == 0)
        model.name = header_value(line);
      else if (line.rfind("norm", 0) == 0)
        model.norm = header_value(line);
      else if (line.rfind("tide_system", 0) == 0)
        model.tide_system = header_value(line);
    }
    if (!found_end || model.GM <= 0.0 || model.R <= 0.0 || model.file_nmax < 0)
      throw std::runtime_error("Incomplete GFC header: " + path);
    if (!model.norm.empty() && model.norm != "fully_normalized")
      throw std::runtime_error("Unsupported normalization: " + model.norm);
    if (nmax < 0)
      nmax = model.file_nmax;
    if (nmax > model.file_nmax)
      throw std::invalid_argument("nmax is larger than the max_degree of the model");
    model.nmax = nmax;
    model.C.assign(tri_size(nmax), 0.0);
    model.S.assign(tri_size(nmax), 0.0);

    // Body: chunks aligned to line starts, one chunk per loop iteration.
    const char *const begin = buf.data() + pos;
    const char *const end = buf.data() + buf.size();
    const std::int64_t body = end - begin;
    const int chunks = static_cast<int>(std::max<std::int64_t>(1, body / (1 << 22)));
    int bad = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : bad)
    for (int ic = 0; ic < chunks; ++ic)
    {
      const char *p = begin + body * ic / chunks;
      const char *stop = begin + body * (ic + 1) / chunks;
      // Move to the first line that starts inside this chunk (chunk 0 starts at a line start).
      if (ic > 0 && p[-1] != '\n')
      {
        p = static_cast<const char *>(std::memchr(p, '\n', end - p));
        p = p ? p + 1 : end;
      }
      while (p < stop)
      {
        const char *eol = static_cast<const char *>(std::memchr(p, '\n', end - p));
        if (!eol)
          eol = end;
        const char *q = p;
        while (q < eol && is_blank(*q))
          ++q;
        // Coefficient line: "gfc" followed by a blank, then n m C S (more fields are ignored).
        if (eol - q > 4 && q[0] == 'g' && q[1] == 'f' && q[2] == 'c' && is_blank(q[3]))
        {
          q += 3;
          int n = 0, m = 0;
          double c = 0, s = 0, sc = 0, ss = 0;
          if (parse_int(q, eol, n) && parse_int(q, eol, m) && parse_real(q, eol, c) && parse_real(q, eol, s) &&
              (!sigma || (parse_real(q, eol, sc) && parse_real(q, eol, ss))) && n >= 0 && m >= 0 && m <= n)
          {
            // Degrees above the requested nmax are parsed and validated but not stored.
            if (n <= nmax)
            {
              const auto i = idx_tri(n, m, nmax);
              model.C[i] = sigma ? sc : c;
              model.S[i] = sigma ? ss : s;
            }
          }
          else
            ++bad;
        }
        p = eol + 1;
      }
    }
    if (bad > 0)
      throw std::runtime_error("GFC with " + std::to_string(bad) + " invalid coefficient lines: " + path);
    return model;
  }
  }
}
