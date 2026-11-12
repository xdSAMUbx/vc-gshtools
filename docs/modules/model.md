# model

Header `src/cpp/gsh_utils/model.hpp`, implementation `src/cpp/gsh_utils/model.cpp`.

The examples in this file use `#include "model.hpp"` (the `gsh` target exposes `src/cpp/gsh_utils` as a public include directory), the `gsh::` namespace, and the EGM2008 file `data/cache/modelos/EGM2008/EGM2008.gfc` relative to the repository root. The reader only builds on Windows as written, so the examples are Windows-only too. The expected outputs are derived by hand from the file header and body, not from a run.

## `struct Model`

| Field | Meaning |
|---|---|
| `GM`, `R` | Gravitational constant (m³/s²) and reference radius (m) |
| `nmax` | Stored degree |
| `file_nmax` | `max_degree` from the header |
| `name`, `norm`, `tide_system` | File metadata |
| `C`, `S` | Coefficients packed with `idx_tri(n, m, nmax)` |

**Example**

```cpp
#include "model.hpp"
#include <iostream>

int main()
{
  std::cout.precision(10);
  gsh::Model empty;                       // defaults: GM = R = 0, nmax = file_nmax = -1
  std::cout << empty.nmax << " " << empty.GM << "\n";   // -1 0

  gsh::Model e = gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
  std::cout << e.name << " " << e.norm << " " << e.tide_system << "\n"; // EGM2008 fully_normalized tide_free
  std::cout << e.GM << " " << e.R << "\n";                  // 3.986004415e+14 6378136.3
  std::cout << e.nmax << " " << e.file_nmax << " " << e.C.size() << "\n"; // 360 2190 65341
}
```

### `double c(int n, int m) const`, `double s(int n, int m) const`
Access to one coefficient. Does not check bounds.

**Example** (`c`)

```cpp
#include "model.hpp"
#include <iostream>

int main()
{
  std::cout.precision(15);
  gsh::Model e = gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
  std::cout << e.c(2, 0) << "\n";  // -0.000484165143790815 (the C20 line of the file)
  std::cout << e.c(1, 0) << "\n";  // 0: EGM2008 has no degree-1 terms, so the slot stays zero
  std::cout << e.c(0, 0) << "\n";  // 1 (stored from "1.0d0")
  // Unchecked: e.c(400, 0) silently returns the slot idx_tri(400, 0, 360), which is the entry for (40, 1).
}
```

**Example** (`s`)

```cpp
#include "model.hpp"
#include <iostream>

int main()
{
  gsh::Model e = gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
  std::cout << e.s(2, 1) << "\n";  // 1.38441389137979e-09 (S21 line of the file)
  std::cout << e.s(2, 0) << "\n";  // 0: S20 is 0.000000000000000e+00 in the file
  std::cout << e.s(0, 0) << "\n";  // 0: stored from "0.0d0"
}
```

### `static Model zeros(double GM, double R, int nmax)`
Model with all coefficients set to 0.

**Example**

```cpp
#include "model.hpp"
#include <iostream>

int main()
{
  // WGS84-like constants: GM = 3.986004418e14 m^3/s^2, R = 6378137 m
  gsh::Model z = gsh::Model::zeros(3.986004418e14, 6378137.0, 10);
  std::cout << z.nmax << " " << z.file_nmax << " " << z.C.size() << "\n"; // 10 10 66 (= tri_size(10))
  std::cout << z.c(10, 10) << " " << z.s(3, 1) << "\n";                    // 0 0
  z.C[gsh::idx_tri(2, 0, 10)] = -4.84165e-4;                                // set C20 by hand
  std::cout << z.c(2, 0) << "\n";                                           // -0.000484165
}
```


### `Model truncated(int new_nmax) const`
Copy limited to a lower degree. Throws `invalid_argument` if `new_nmax > nmax`.

**Example**

```cpp
#include "model.hpp"
#include <iostream>
#include <stdexcept>

int main()
{
  gsh::Model full = gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
  gsh::Model low = full.truncated(10);
  std::cout << low.nmax << " " << low.file_nmax << " " << low.C.size() << "\n"; // 10 2190 66
  std::cout << (low.c(2, 0) == full.c(2, 0)) << "\n";                          // 1 (same coefficient)
  try
  {
    gsh::Model bad = full.truncated(400);   // 400 > 360
  }
  catch (const std::invalid_argument &e)
  {
    std::cout << e.what() << "\n"; // truncated: degree larger than the model degree
  }
}
```

## `Model load_gfc(const std::string &path, int nmax = -1)`

Reads an ICGEM `.gfc` file. With EGM2008 it takes 0.10 s.

- **`nmax`:** with `nmax < 0` uses `max_degree`; a larger value throws `invalid_argument`.
- **Header:** read in sequence up to `end_of_head`. Requires `earth_gravity_constant`, `radius` and `max_degree`.
- **Normalization:** a `norm` other than `fully_normalized` throws `runtime_error` ("Unsupported normalization"). A file with no `norm` line is accepted and treated as fully normalized.
- **Body:** split into chunks of about 4 MiB aligned to line starts, and each thread processes one chunk.
- **Missing coefficients:** remain at 0. For example, EGM2008 does not include degree 1.
- **Fortran exponents:** accepts `1.0d0`.
- **Errors:** throws `runtime_error` if the file cannot be opened or read completely, the header is incomplete or has a malformed numeric value, the normalization is unsupported, or there are invalid `gfc` lines (it reports how many).

As written, it only builds on Windows, because it uses `_fseeki64` and `_ftelli64`.

**Example**

```cpp
#include "model.hpp"
#include <iostream>

int main()
{
  // EGM2008 declares max_degree 2190, so nmax = -1 keeps every degree up to 2190.
  gsh::Model egm = gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc");
  std::cout << egm.nmax << " " << egm.C.size() << "\n";  // 2190 2401336 (= tri_size(2190), about 19 MB per array)
  std::cout << egm.c(2, 0) << "\n";                     // -0.000484165143790815
  std::cout << egm.c(0, 0) << "\n";                     // 1
}
```

**Example** (reduced degree)

```cpp
#include "model.hpp"
#include <iostream>

int main()
{
  // Read only degrees <= 360: far fewer coefficients to store and fewer to parse.
  gsh::Model m = gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
  std::cout << m.nmax << " " << m.file_nmax << " " << m.C.size() << "\n"; // 360 2190 65341
}
```

**Example** (errors)

```cpp
#include "model.hpp"
#include <iostream>
#include <stdexcept>

int main()
{
  try
  {
    gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 3000); // 3000 > max_degree 2190
  }
  catch (const std::invalid_argument &e)
  {
    std::cout << "invalid_argument: " << e.what() << "\n"; // nmax is larger than the max_degree of the model
  }
  try
  {
    gsh::load_gfc("data/cache/modelos/EGM2008/missing.gfc"); // no such file
  }
  catch (const std::runtime_error &e)
  {
    std::cout << "runtime_error: " << e.what() << "\n"; // Cannot open: data/cache/modelos/EGM2008/missing.gfc
  }
  // A file with a malformed "gfc" line throws runtime_error ("GFC with N invalid coefficient lines: ...").
  // The repo has no such file, so that path is not run here.
}
```

## Internal functions (`model.cpp`)

| Function | What it does |
|---|---|
| `read_file(path)` | Loads the whole file into memory |
| `is_blank(c)` | True for space, tab or `\r` |
| `parse_real(p, end, out)` | Reads a real and advances `p`. Tries `from_chars`; if that fails, replaces `d/D` with `e` and retries |
| `parse_int(p, end, out)` | Reads an integer and advances `p` |
| `header_value(line)` | Returns the second token of a header line |

These helpers are in an anonymous namespace in `model.cpp`, so user code cannot call them. Each example below shows the behavior through `load_gfc`, and each one says what it needs.

**`read_file(path)`**

**Example**

```cpp
#include "model.hpp"
#include <iostream>
#include <stdexcept>

int main()
{
  // read_file throws here because fopen fails. The "Incomplete read" branch needs a truncated file, which the repo does not have.
  try
  {
    gsh::load_gfc("data/cache/modelos/EGM2008/no_such_file.gfc");
  }
  catch (const std::runtime_error &e)
  {
    std::cout << e.what() << "\n"; // Cannot open: data/cache/modelos/EGM2008/no_such_file.gfc
  }
}
```

**`is_blank(c)`**

**Example**

```cpp
#include "model.hpp"
#include <iostream>

int main()
{
  // is_blank treats ' ', '\t' and '\r' as separators inside a data line, so a CRLF-saved .gfc loads.
  // The branch needs a copy of the file with CRLF endings or tab separators, which is not in the repo.
  // The call is the same for any such copy:
  gsh::Model m = gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
  std::cout << m.c(2, 0) << "\n"; // -0.000484165143790815 for any line ending
}
```

**`parse_real(p, end, out)`**

**Example**

```cpp
#include "model.hpp"
#include <iostream>

int main()
{
  // The degree-0 line is "gfc 0 0 1.0d0 0.0d0 ...". from_chars stops at 'd', so parse_real retries
  // with 'd' replaced by 'e' and stores 1.0.
  gsh::Model m = gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 2);
  std::cout << m.c(0, 0) << " " << m.s(0, 0) << "\n"; // 1 0
}
```

**`parse_int(p, end, out)`**

**Example**

```cpp
#include "model.hpp"
#include <iostream>

int main()
{
  // parse_int reads the n and m fields of "gfc 2 1 -0.206615509074176e-09 0.138441389137979e-08 ...".
  // Those values place the coefficients at idx_tri(2, 1, nmax).
  gsh::Model m = gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 2);
  std::cout << m.c(2, 1) << " " << m.s(2, 1) << "\n"; // -2.06615509074176e-10 1.38441389137979e-09
}
```

**`header_value(line)`**

**Example**

```cpp
#include "model.hpp"
#include <iostream>

int main()
{
  // header_value returns the second token: "modelname    EGM2008" gives "EGM2008".
  gsh::Model m = gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 2);
  std::cout << m.name << " | " << m.norm << " | " << m.tide_system << "\n";
  // EGM2008 | fully_normalized | tide_free
}
```

## Putting it together

A typical workflow: load a model at a chosen degree, check its metadata, reduce it for a quick check, build a synthetic model with the same constants, and handle the errors that the reader throws.

```cpp
#include "model.hpp"
#include <iostream>
#include <stdexcept>

int main()
{
  std::cout.precision(10);
  try
  {
    // 1. Load EGM2008 up to degree 360 (about 65k coefficients per array instead of 2.4 M).
    gsh::Model full = gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 360);
    std::cout << full.name << ", GM = " << full.GM << " m^3/s^2, R = " << full.R << " m\n";
    std::cout << "stored degree " << full.nmax << " of file degree " << full.file_nmax << "\n";

    // 2. Keep only degrees <= 20 for a faster run.
    gsh::Model low = full.truncated(20);
    std::cout << "C20 = " << low.c(2, 0) << "\n"; // -0.0004841651438

    // 3. Build a synthetic model with only C20 set, using the same GM and R.
    gsh::Model test = gsh::Model::zeros(full.GM, full.R, 20);
    test.C[gsh::idx_tri(2, 0, 20)] = low.c(2, 0);
    std::cout << "synthetic C20 = " << test.c(2, 0) << ", C30 = " << test.c(3, 0) << "\n"; // -0.0004841651438, 0

    // 4. A degree beyond the file's max_degree is rejected.
    gsh::Model bad = gsh::load_gfc("data/cache/modelos/EGM2008/EGM2008.gfc", 5000);
  }
  catch (const std::invalid_argument &e)
  {
    std::cerr << "bad degree: " << e.what() << "\n";
  }
  catch (const std::runtime_error &e)
  {
    std::cerr << "file problem: " << e.what() << "\n";
  }
}
```
