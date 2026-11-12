# table_io

Header `src/cpp/gsh_utils/table_io.hpp`, implementation `src/cpp/gsh_utils/table_io.cpp`.

## `Table read_table(path, sep, header_lines, cols)`

Reads numeric columns from a text file.

| Parameter | Meaning |
|---|---|
| `path` | File |
| `sep` | Delimiter. `'\0'` means runs of spaces or tabs |
| `header_lines` | Leading lines to ignore |
| `cols` | Column indices, zero-based. The result follows this order |

**Returns** `Table{columns, skipped_lines}`. `columns[c][i]` is the value of the requested column `c` in valid row `i`. `skipped_lines` lists the skipped line numbers, one-based.

**Invalid rows:** if a requested field is not numeric or fields are missing, the whole row is skipped and recorded. It is never corrected or interpolated. This is how line 13 of `GRAVD_ms05_median_debiased_1hz.txt` is handled. Blank lines are ignored without being recorded.

**Errors:** throws `runtime_error` if the file cannot be opened or the read is incomplete.

All examples below assume `#include "table_io.hpp"` (with `src/cpp/gsh_utils` on the include path, since the header includes `core.hpp` from the same folder), the `gsh` namespace, and that the program runs from the repository root so the relative data paths resolve.

**Example** (basic read, Colorado-style file):

```cpp
#include <cstdio>
#include "table_io.hpp"

int main()
{
  using namespace gsh;
  // Extracted copy of the Colorado gravity file. Fields: block,line,time,lat,lon,h,g.
  // Skip the header line and keep lat (3), lon (4), h (5), g (6), in that order.
  Table t = read_table("data/cache/tablas/GRAVD_ms05_median_debiased_1hz/GRAVD_ms05_median_debiased_1hz.txt",
                       ',', 1, {3, 4, 5, 6});
  std::printf("valid rows=%zu skipped=%zu\n", t.columns[0].size(), t.skipped_lines.size());
  // First valid row: lat=38.20304581 deg, lon=-105.65179298 deg, h=5376.261 m, g=978431.7755 mGal
  std::printf("first: %.8f %.8f %.3f %.4f\n",
              t.columns[0][0], t.columns[1][0], t.columns[2][0], t.columns[3][0]);
}
```

**Example** (invalid rows are skipped, not repaired). Sample file `bad.csv`, which is not part of the repository:

```text
x,y,z
1.0,2.0,3.0
4.0,oops,6.0
7.0,8.0
+9.5,10,11
```

```cpp
#include <cstdio>
#include "table_io.hpp"

int main()
{
  using namespace gsh;
  // Line 1 is the header (header_lines = 1). Requested fields: x, y, z.
  Table t = read_table("bad.csv", ',', 1, {0, 1, 2});
  // Line 3 ("oops" is not a number) and line 4 (field z is missing) are rejected.
  // Line 5 is accepted; the leading '+' is allowed.
  for (std::int64_t ln : t.skipped_lines)
    std::printf("rejected line %lld\n", static_cast<long long>(ln)); // rejected line 3, rejected line 4
  // Only requested fields are checked: with cols {0, 2}, line 3 would be kept (y is not read).
}
```

**Example** (the `Table` result). Uses `bad.csv` from above:

```cpp
#include <cstdio>
#include "table_io.hpp"

int main()
{
  using namespace gsh;
  Table t = read_table("bad.csv", ',', 1, {0, 1, 2});
  // t.columns has one vector per requested column, all with the same length (here 2).
  for (std::size_t i = 0; i < t.columns[0].size(); ++i)
    std::printf("x=%.1f y=%.1f z=%.1f\n", t.columns[0][i], t.columns[1][i], t.columns[2][i]);
  // Output: x=1.0 y=2.0 z=3.0, then x=9.5 y=10.0 z=11.0
}
```

**Example** (error handling). The file does not exist, so `read_table` throws:

```cpp
#include <cstdio>
#include <stdexcept>
#include "table_io.hpp"

int main()
{
  using namespace gsh;
  try
  {
    read_table("data/does_not_exist.txt", ',', 0, {0});
  }
  catch (const std::runtime_error &e)
  {
    std::fprintf(stderr, "%s\n", e.what()); // prints: Cannot open: data/does_not_exist.txt
    return 1;
  }
  return 0;
}
```

## Use with the repository data

```cpp
// Colorado: lat, lon, h, g
read_table("data/cache/tablas/GRAVD_ms05_median_debiased_1hz/GRAVD_ms05_median_debiased_1hz.txt", ',', 1, {3, 4, 5, 6});
// Italy: Lon, Lat, h, perturbation. Longitude comes first.
read_table("data/italia/datos_aereos/airborne_data.txt", '\0', 1, {0, 1, 2, 3});
```

**Note on the file formats.** The Colorado file is distributed zipped in `data/colorado/` (`GRAVD_ms05_median_debiased_1hz.zip`) and the examples use the extracted copy under `data/cache/tablas/`. The Colorado file is comma-separated with one header line:

```text
block,line,time(YYYYMMDDHHMMSS.SSS),lat(deg),lon(deg),h(m),g (mGal)
MS05,103,20170823150118400,38.20304581,-105.65179298,5376.261,978431.7755
MS05,103,20170823150119400,38.20307974,-105.65081569,5375.812,978431.229
```

The Italian airborne file is whitespace-separated (spaces or tabs) with one header line:

```text
Lon	Lat	h	gravity disturbance(mGal)
15.308	44.992	2500	52.50767444878
15.377	44.948	2500	58.95392153068
```

The extracted Colorado copy parses line 13 without problems, so the statement above about line 13 should be checked against the original zip before it is relied on.

**Example** (Italian airborne file, whitespace mode):

```cpp
#include <cstdio>
#include "table_io.hpp"

int main()
{
  using namespace gsh;
  Table it = read_table("data/italia/datos_aereos/airborne_data.txt", '\0', 1, {0, 1, 2, 3});
  // Row 0: Lon=15.308 deg, Lat=44.992 deg, h=2500 m, gravity disturbance=52.50767444878 mGal
  std::printf("rows=%zu skipped=%zu\n", it.columns[0].size(), it.skipped_lines.size());
  std::printf("first: lon=%.3f lat=%.3f h=%.0f dg=%.5f\n",
              it.columns[0][0], it.columns[1][0], it.columns[2][0], it.columns[3][0]);
}
```

## Putting it together

A typical workflow: read the Colorado table, report the rejected lines, and summarise the ellipsoidal heights. Only `read_table` and the standard library are used.

```cpp
#include <algorithm>
#include <cstdio>
#include <numeric>
#include <stdexcept>
#include "table_io.hpp"

int main()
{
  using namespace gsh;
  try
  {
    // Keep lat (3), lon (4), h (5), g (6) from the extracted Colorado copy.
    const std::string path =
        "data/cache/tablas/GRAVD_ms05_median_debiased_1hz/GRAVD_ms05_median_debiased_1hz.txt";
    Table t = read_table(path, ',', 1, {3, 4, 5, 6});

    // Report what was rejected, using 1-based line numbers as shown in an editor.
    std::printf("valid=%zu skipped=%zu\n", t.columns[0].size(), t.skipped_lines.size());
    for (std::int64_t ln : t.skipped_lines)
      std::printf("  rejected line %lld\n", static_cast<long long>(ln));

    // Summary of the ellipsoidal height h (column 2, metres).
    const std::vector<double> &h = t.columns[2];
    if (!h.empty())
    {
      const double mean = std::accumulate(h.begin(), h.end(), 0.0) / static_cast<double>(h.size());
      const auto [mn, mx] = std::minmax_element(h.begin(), h.end());
      std::printf("h: mean=%.3f min=%.3f max=%.3f m\n", mean, *mn, *mx);
    }
  }
  catch (const std::runtime_error &e)
  {
    std::fprintf(stderr, "read failed: %s\n", e.what());
    return 1;
  }
  return 0;
}
```
