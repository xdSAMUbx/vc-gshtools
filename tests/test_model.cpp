// test_model.cpp: tests of the coefficient model reader (load_gfc, Model::truncated) on a small
// synthetic .gfc file and on the EGM2008 file when it is present, plus the delimited table reader
// (read_table) on a synthetic file with one corrupt line.
#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include "model.hpp"
#include "table_io.hpp"
using namespace gsh;

namespace fs = std::filesystem;

// Small synthetic file: header fields, Fortran exponents (1.0d0, 2.5D-6) and extra trailing fields.
TEST(Model, ParsesGfcSmall)
{
  const auto path = fs::temp_directory_path() / "gsh_test.gfc";
  {
    std::ofstream o(path);
    o << "modelname X\nearth_gravity_constant 0.3986004415E+15\nradius 0.63781363E+07\nmax_degree 3\n"
         "norm fully_normalized\ntide_system tide_free\nend_of_head ====\n"
         "gfc 0 0 1.0d0 0.0d0 0 0\ngfc 2 0 -1.5e-3 0 0 0\ngfc 3 2 2.5D-6 -3.5e-7 0 0\n";
  }
  const auto m = load_gfc(path.string());
  EXPECT_EQ(m.nmax, 3);
  EXPECT_DOUBLE_EQ(m.GM, 0.3986004415e15);
  EXPECT_DOUBLE_EQ(m.R, 6378136.3);
  EXPECT_DOUBLE_EQ(m.c(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(m.c(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(m.c(2, 0), -1.5e-3);
  EXPECT_DOUBLE_EQ(m.c(3, 2), 2.5e-6);
  EXPECT_DOUBLE_EQ(m.s(3, 2), -3.5e-7);
  const auto t = m.truncated(2);
  EXPECT_EQ(t.nmax, 2);
  EXPECT_DOUBLE_EQ(t.c(2, 0), -1.5e-3);
  EXPECT_THROW(load_gfc(path.string(), 9), std::invalid_argument);
  fs::remove(path);
}

// Reads the EGM2008 file from the repository data folder if it exists; skipped otherwise.
TEST(Model, ReadsEgm2008IfPresent)
{
  const std::string path = std::string(GSH_REPO_DATA) + "/cache/modelos/EGM2008/EGM2008.gfc";
  if (!fs::exists(path))
    GTEST_SKIP();
  const auto m = load_gfc(path, 10);
  EXPECT_DOUBLE_EQ(m.c(2, 0), -0.484165143790815e-03);
  EXPECT_DOUBLE_EQ(m.s(5, 5), -0.669379935180165e-06);
  EXPECT_EQ(m.tide_system, "tide_free");
}

// A non-numeric field rejects its whole line (reported in skipped_lines, not repaired).
TEST(TableIo, SkipsCorruptLineWithoutGuessing)
{
  const auto path = fs::temp_directory_path() / "gsh_table.txt";
  {
    std::ofstream o(path);
    o << "a,b,c\nX,1.5,2\nX,<,3\nX,2.5,4\n";
  }
  const auto t = read_table(path.string(), ',', 1, {1, 2});
  ASSERT_EQ(t.columns[0].size(), 2u);
  EXPECT_DOUBLE_EQ(t.columns[0][1], 2.5);
  ASSERT_EQ(t.skipped_lines.size(), 1u);
  EXPECT_EQ(t.skipped_lines[0], 3);
  fs::remove(path);
}
