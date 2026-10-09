#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
build_dir="${SHIP_BUILD_DIR:-${repo_dir}/build-ship}"
cmake_bin="${CMAKE_BIN:-cmake}"
"${cmake_bin}" -S "${repo_dir}" -B "${build_dir}" \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_UNITTESTS=OFF -DBUILD_BENCHMARKS=OFF \
  -DBUILD_EXAMPLES=ON -DWITH_OPENMP=OFF -DGIT_SUBMOD_AUTO=OFF \
  -DNATIVE_SIZE=64 -DCMAKE_POLICY_VERSION_MINIMUM=3.5
"${cmake_bin}" --build "${build_dir}" --parallel "${SHIP_JOBS:-4}" \
  --target ship-baseline-checks ship-component-checks ship-half-bootstrap-checks
for name in ship-baseline-checks ship-component-checks ship-half-bootstrap-checks; do
  "${build_dir}/bin/examples/pke/${name}"
done
