# Lattice-estimator check of the SHIP parameter sets.
#
# Run with SageMath and https://github.com/malb/lattice-estimator on the path, e.g.
#   sage -python research/ship/security/ship-estimator.sage
# The results in research/ship/security/estimator-results.txt were produced with
# lattice-estimator commit d8c00b483f70b6298690386b1a3278b5167a2e42 (2026-08-19) and SageMath 10.5
# (mybinder.org image of malb/lattice-estimator).
#
# Dense keys: OpenFHE uniform ternary secret, every evaluation key lives modulo Q*P.
# Sparse key: h = 31 ternary secret, used only by the dense-to-sparse switching key modulo q0 * p'.
# The sparse distribution is modelled as SparseTernary(16, 15, n) (16 coefficients +1, 15 coefficients -1).
# The position restriction of SHIP (one nonzero coefficient per window of width 2w) is not modelled by the
# estimator; it is covered separately by the [May21] bound printed by ship-paper-bench.

from estimator import *
import math
import time

Logging.set_level(Logging.LEVEL0)


def log2rop(result):
    return float(math.log2(result["rop"]))


def run(tag, n, logq, sparse, m=oo):
    xs = ND.SparseTernary(16, 15, n) if sparse else ND.UniformMod(3, n)
    params = LWE.Parameters(n=n, q=2**logq, Xs=xs, Xe=ND.DiscreteGaussian(3.19), m=m)
    attacks = [("usvp", LWE.primal_usvp), ("dual_hybrid", LWE.dual_hybrid)]
    if sparse:
        # The attack LWE.estimate() reports as bdd_mitm_hybrid.
        attacks.append(("bdd_mitm_hybrid", lambda p: LWE.primal_hybrid(p, mitm=True, babai=True)))
    for name, attack in attacks:
        start = time.time()
        try:
            value = "%.1f" % log2rop(attack(params))
        except Exception as e:  # noqa: BLE001
            value = "error: %s" % e
        print("%-6s %-6s N=2^%d log2q=%d m=%s | %-16s | %s bits | %.0f s"
              % (tag, "sparse" if sparse else "dense", int(math.log2(n)), logq, "n" if m == n else "inf",
                 name, value, time.time() - start), flush=True)


# Dense keys (log2(QP) of the contexts) and the HE-standard reference points.
for tag, n, logq in [("ref", 8192, 218), ("LL13", 8192, 218), ("ref", 16384, 438), ("LL14", 16384, 437),
                     ("HT14", 16384, 434), ("ref", 32768, 881), ("HT15", 32768, 827)]:
    run(tag, n, logq, sparse=False)

# Sparse key (log2(q0 p') bounds used by the presets).
for tag, n, logq in [("LL13", 8192, 55), ("HT14", 16384, 88), ("LL14", 16384, 100), ("HT15", 32768, 105)]:
    run(tag, n, logq, sparse=True)
    run(tag, n, logq, sparse=True, m=n)
