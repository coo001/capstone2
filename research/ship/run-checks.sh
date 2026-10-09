#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CMAKE_BIN="${CMAKE_BIN:-cmake}"
SHIP_BUILD_DIR="${SHIP_BUILD_DIR:-${ROOT}/build-ship}"
SHIP_JOBS="${SHIP_JOBS:-4}"

"${CMAKE_BIN}" -S "${ROOT}" -B "${SHIP_BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_BENCHMARKS=OFF -DBUILD_UNITTESTS=OFF \
    -DBUILD_EXAMPLES=ON -DWITH_OPENMP=OFF -DGIT_SUBMOD_AUTO=OFF \
    -DNATIVE_SIZE=64 -DCMAKE_POLICY_VERSION_MINIMUM=3.5

targets=(ship-baseline-checks ship-component-checks ship-fused-rotation-checks ship-half-bootstrap-checks ship-full-bootstrap-checks)
"${CMAKE_BIN}" --build "${SHIP_BUILD_DIR}" --target "${targets[@]}" --parallel "${SHIP_JOBS}"
for target in "${targets[@]}"; do
    "${SHIP_BUILD_DIR}/bin/examples/pke/${target}"
done
