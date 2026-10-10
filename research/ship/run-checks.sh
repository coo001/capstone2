#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CMAKE_BIN="${CMAKE_BIN:-cmake}"
SHIP_BUILD_DIR="${SHIP_BUILD_DIR:-${ROOT}/build-ship}"
SHIP_JOBS="${SHIP_JOBS:-4}"

# WITH_REDUCED_NOISE enables centered rounding in hybrid key switching (required for the
# precision of the 128-bit SHIP parameter sets; see research/ship/LIBRARY.md).
"${CMAKE_BIN}" -S "${ROOT}" -B "${SHIP_BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_BENCHMARKS=OFF -DBUILD_UNITTESTS=ON \
    -DBUILD_EXAMPLES=ON -DWITH_OPENMP=OFF -DWITH_REDUCED_NOISE=ON -DGIT_SUBMOD_AUTO=OFF \
    -DNATIVE_SIZE=64 -DCMAKE_POLICY_VERSION_MINIMUM=3.5

# Research prototypes (examples/pke/ship/*.h) and the library implementation (lib/scheme/ckksrns/ckksrns-ship.cpp).
targets=(ship-baseline-checks ship-component-checks ship-fused-rotation-checks ship-half-bootstrap-checks
         ship-full-bootstrap-checks ship-aux-masking-checks ship-library-checks pke_tests)
"${CMAKE_BIN}" --build "${SHIP_BUILD_DIR}" --target "${targets[@]}" --parallel "${SHIP_JOBS}"
for target in "${targets[@]}"; do
    if [ "${target}" = pke_tests ]; then
        "${SHIP_BUILD_DIR}/unittest/pke_tests" --gtest_filter='UTCKKSRNS_SHIP*'
    else
        "${SHIP_BUILD_DIR}/bin/examples/pke/${target}"
    fi
done
# 128-bit parameter sets (about 10 GB / 16 GB of memory):
#   "${SHIP_BUILD_DIR}/bin/examples/pke/ship-paper-bench" LL13 5
#   "${SHIP_BUILD_DIR}/bin/examples/pke/ship-paper-bench" LL14 5
