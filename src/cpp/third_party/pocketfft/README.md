# pocketfft (vendored)

pocketfft is a header-only C++ FFT library. It is vendored here without modification and is used only by
`src/cpp/gsh_expansion/grid.cpp` for the real backward transform of the longitude sums on full-circle grids.

- Upstream: https://github.com/mreineck/pocketfft, branch `cpp`
- Commit taken: `4d865cba24359b05b37da4e75aee24498c3b78e7` (committed 2026-10-18, "Add POCKETFFT_NO_EXCEPTIONS option
  for builds without exception handling"). The header has no version number; the commit identifies it.
- Files: `pocketfft_hdronly.h` (the upstream name; the brief assumed `.hpp`), `LICENSE.md` (verbatim upstream copy).
- SHA-256 of `pocketfft_hdronly.h`: `bc0095b216ec71fd66686e6b7f4d572c2d1af7f2f44848dde8071426f5d5f16d`.
- Licence: BSD-3-Clause, see `LICENSE.md`. Copyright (C) 2010-2018 Max-Planck-Society; authors Martin Reinecke and Peter Bell.

Build notes: `grid.cpp` defines `POCKETFFT_NO_MULTITHREADING` before the include, so the header does not use
`std::thread` or pthreads; the grid is parallelised with OpenMP instead. The plan (`pocketfft_r<double>`) is built
once per synthesis and its `exec()` is called from several threads (it keeps its scratch on the call stack).

To update: replace the header and `LICENSE.md` from upstream, record the new commit here, and rerun the grid tests
(`Grid.C2RMatchesNaiveDft` in `tests/cpp/test_grid.cpp` checks the transform convention).
