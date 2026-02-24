# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Earcut.hpp is a header-only C++14 polygon triangulation library — a C++ port of Mapbox's earcut.js. The entire library lives in a single header: `include/mapbox/earcut.hpp`. It uses a modified ear slicing algorithm with z-order curve hashing for performance.

## Build Commands

```bash
# Configure (from repo root, creates build/ directory)
mkdir -p build && cd build && cmake .. -GNinja

# Build all targets (tests, benchmarks, visualizer)
ninja                          # from build/

# Run tests
build/tests                    # from repo root, or ./tests from build/

# Run benchmarks
build/bench

# Build without optional targets
cmake .. -GNinja -DEARCUT_BUILD_VIZ=OFF -DEARCUT_BUILD_BENCH=OFF

# Format check (CI uses clang-format)
find include test -name "*.hpp" -o -name "*.cpp" | xargs clang-format --dry-run --Werror

# Auto-format
find include test -name "*.hpp" -o -name "*.cpp" | xargs clang-format -i
```

Tests use GoogleTest with parametrized fixture tests. There is no way to run a single fixture in isolation without using `--gtest_filter`:
```bash
build/tests --gtest_filter="FixtureTests/EarcutAreaTest.EarcutTriangulation/water"
```

## Architecture

**Library** (`include/mapbox/earcut.hpp`): Single-header, template-based. The public API is:
- `mapbox::earcut<N>(polygon)` — free function returning `std::vector<N>` of triangle indices
- `mapbox::Earcut<N>` — class with reusable `operator()`, stores results in `indices` and `vertices` members
- `mapbox::util::nth<I, T>` — point coordinate accessor, specialized for `std::tuple`, `std::pair`, `std::array`; extend this for custom point types

Key internals: `Node` (doubly-linked list vertex with z-order pointers), `ObjectPool<T>` (block-based memory pool reusable across runs), `Triangle` (cache-optimized containment test struct).

**Tests** (`test/`):
- `test/test.cpp` — GoogleTest parametrized tests comparing earcut output (triangle count and area deviation) against expected values and libtess2
- `test/bench.cpp` — Google Benchmark suite over selected fixture geometries
- `test/viz.cpp` — OpenGL interactive visualizer (requires GLFW3 + OpenGL)
- `test/fixtures/*.cpp` — Each file defines a polygon geometry with expected triangle count and deviation tolerances, auto-registered via `Collector<T>` singleton
- `test/fixtures/geometries.hpp` — `FixtureTester`/`Fixture<T>` framework for test parametrization
- `test/comparison/` — Wrapper interfaces for earcut and libtess2 used by the test framework

**Vendor** (`vendor/`): Git submodules — googletest, benchmark, glfw. Initialize with `git submodule update --init --recursive`.

## Adding a Test Fixture

Create a new `.cpp` file in `test/fixtures/` following the pattern of existing fixtures: define a `Fixture<CoordType>` with name, expected triangle count, earcut deviation tolerance, libtess2 deviation tolerance, and polygon data. The fixture auto-registers itself and will be picked up by the parametrized test suite.

## Key Build Details

- C++14 required
- Floating-point precision: fixtures and tests compile with `-ffp-contract=off` (non-MSVC) to ensure consistent cross-platform results
- UBSan enabled automatically in Debug builds (Clang/GCC)
- CI runs: Linux (gcc/clang), macOS, sanitizers (ASan+UBSan), coverage (lcov), and clang-format checks
- Code style: `.clang-format` based on Google style, 120 column limit, 4-space indent
