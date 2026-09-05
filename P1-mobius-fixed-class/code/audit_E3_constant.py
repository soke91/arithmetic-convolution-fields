# -*- coding: utf-8 -*-
r"""
paper/theorem_A.md, Section {#sec:C}, item (i) -- the sign of the residual
main term, and of the constant that fixes it.

WHAT IS UNDER TEST

The section evaluates the log-weighted functional

    E_3(alpha) = sum_{k<K,(k,N)=1} mu(k) log k * E_mu(N;k)

by the split T_w = P_w - MT_w - C(N)*B_w of [eq:split0]/[eq:R3], and prints

    (i)  MT_log = -S(N)*N,   "the constant is fixed by A(N)*Gt(1) = S(N)"
    (ii) B_log(K) -> -S(N)

with Gt(1) = lim_x sum_{m<=x} mu(m) lambda(m) 1_{(m,N)=1} log m / m,
lambda(m) = prod_{p|m} (1 - 1/(p(p-1)))^{-1},
A(N) = prod_{p not| N, p>2} (1 - 1/(p(p-1))),
S(N) = 2 prod_{p>2}(1-1/(p-1)^2) prod_{p|N,p>2}(1+1/(p-2)).

Chaining (i) and (ii) through T = P - MT - C*B gives

    E_3 = rt(N) - MT_log + S(N) C(N),

so the printed MT_log = -S(N)N yields E_3 = rt(N) + S(N)N + S(N)C(N),
while Theorem [thm:C] on the same page asserts

    E_3 = rt(N) - S(N)(N - C(N)) + O_A(N (log N)^{-A}).

The two differ by 2*S(N)*N, which is of the size of the whole object.
At most one of them is right.

PRE-REGISTERED PREDICTION (written before this script was run)

  P1  A(N)*Gt(1) converges to -S(N), not to +S(N).  Equivalently the
      printed constant identity carries the wrong sign and MT_log = +S(N)N.

  P2  B_log(K) converges to -S(N).  (The paper's (ii); expected to hold.)

  P3  The brute-force E_3 at finite N sits near rt(N) - S(N)(N - C(N))
      and nowhere near rt(N) + S(N)N + S(N)C(N).

REFUTATION RULE (fixed before the run)

  P1 is REFUTED if A(N)*Gt(1) is positive at every x tested, or if
     |A(N)Gt(1) - S(N)| < |A(N)Gt(1) + S(N)| at the largest x.
  P2 is REFUTED if B_log(K) is positive at every K tested.
  P3 is REFUTED if |E_3 - (rt - S(N-C))| > |E_3 - (rt + S*N + S*C)| at
     the largest N tested.
  The script exits non-zero if any of P1, P2, P3 is refuted -- which is
  the outcome that says the paper is right and this audit is wrong.

WHY THE FINITE-N TEST IS DECISIVE DESPITE A LARGE ERROR TERM

Theorem [thm:C]'s error is O_A(N(log N)^{-A}), which at N = 10^6 is not
small. But the two candidate right-hand sides differ by 2*S(N)*N ~ 3.5e6,
an order of magnitude above anything the error term can hide at these N.
The test separates the sign; it does not test the error term.

CITED BY: {#rem:sign} in paper/.
"""

import io
import math
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "audit_E3_constant.txt")

THETA = 0.56
TWIN = None            # 2 * prod_{p>2} (1 - 1/(p-1)^2), computed below


def sieve(n):
    """mu[0..n], Lambda[0..n], smallest-prime-factor based."""
    spf = np.zeros(n + 1, dtype=np.int64)
    for p in range(2, n + 1):
        if spf[p] == 0:
            spf[p::p] = np.where(spf[p::p] == 0, p, spf[p::p])
    mu = np.ones(n + 1, dtype=np.int64)
    mu[0] = 0
    lam = np.zeros(n + 1, dtype=np.float64)
    for v in range(2, n + 1):
        p = spf[v]
        w = v // p
        if w % p == 0:
            mu[v] = 0
        else:
            mu[v] = -mu[w]
        lam[v] = math.log(p) if (w == 1 or spf[w] == p and _ppow(v, p)) else 0.0
    return mu, lam, spf


def _ppow(v, p):
    while v % p == 0:
        v //= p
    return v == 1


def von_mangoldt(n, spf):
    lam = np.zeros(n + 1, dtype=np.float64)
    for p in range(2, n + 1):
        if spf[p] != p:
            continue
        q = p
        while q <= n:
            lam[q] = math.log(p)
            if q > n // p:
                break
            q *= p
    return lam


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.flatnonzero(s)


def twin_constant(limit=4_000_000):
    """2 * prod_{p>2} (1 - 1/(p-1)^2), with the tail estimated."""
    ps = primes_upto(limit)
    val = 2.0
    for p in ps[1:]:
        val *= 1.0 - 1.0 / (p - 1.0) ** 2
    # tail: log-prod over p>limit is ~ -sum 1/(p-1)^2 ~ -1/(limit log limit)
    return val


def singular_series(N, ps_of_N):
    v = TWIN
    for p in ps_of_N:
        if p > 2:
            v *= 1.0 + 1.0 / (p - 2.0)
    return v


def A_of_N(N, ps_of_N, limit=4_000_000):
    """prod_{p not| N, p>2} (1 - 1/(p(p-1)))."""
    ps = primes_upto(limit)
    bad = set(int(p) for p in ps_of_N)
    v = 1.0
    for p in ps[1:]:
        if int(p) in bad:
            continue
        v *= 1.0 - 1.0 / (p * (p - 1.0))
    return v


def prime_factors(N, spf):
    out = []
    while N > 1:
        p = int(spf[N])
        out.append(p)
        while N % p == 0:
            N //= p
    return out


def G_tilde(x, mu, spf, ps_of_N):
    """sum_{m<=x} mu(m) lambda(m) 1_{(m,N)=1} log m / m, partial sums."""
    bad = set(ps_of_N)
    tot = 0.0
    marks = {}
    checkpoints = [x // 64, x // 16, x // 4, x]
    for m in range(2, x + 1):
        # The checkpoint is a place in the sum, not a term of it, so it
        # is recorded BEFORE the two skips.  Recorded after them, an
        # index that is not squarefree or shares a prime with N is
        # never written: at N = 10^6 and 4*10^6 all three of x//64,
        # x//16, x//4 carry 4 or 25, so marks held one point and the
        # registered clause "positive at every x tested" was scored on
        # that one point.
        if m in checkpoints:
            marks[m] = tot
        if mu[m] == 0:
            continue
        v, lam, ok = m, 1.0, True
        while v > 1:
            p = int(spf[v])
            if p in bad:
                ok = False
                break
            lam /= (1.0 - 1.0 / (p * (p - 1.0)))
            while v % p == 0:
                v //= p
        if not ok:
            continue
        tot += mu[m] * lam * math.log(m) / m
    marks[x] = tot
    return marks


def B_log(K, mu, spf, ps_of_N):
    """sum_{k<K,(k,N)=1} mu(k) log k / phi(k), partial sums."""
    bad = set(ps_of_N)
    tot = 0.0
    marks = {}
    checkpoints = [K // 64, K // 16, K // 4, K - 1]
    for k in range(2, K):
        # Recorded before the skips, for the reason given in G_tilde.
        # Here the field happened to leave two or three checkpoints
        # squarefree and coprime, so the table was not empty -- the
        # hole was the same one and only the field differed.
        if k in checkpoints:
            marks[k] = tot
        if mu[k] == 0:
            continue
        v, phi, ok = k, 1, True
        while v > 1:
            p = int(spf[v])
            if p in bad:
                ok = False
                break
            phi *= (p - 1)
            while v % p == 0:
                v //= p
        if not ok:
            continue
        tot += mu[k] * math.log(k) / phi
    marks[K - 1] = tot
    return marks


def E3_bruteforce(N, mu, lam, spf, K):
    """E_3(alpha) = sum_{k<K,(k,N)=1} mu(k) log k * E_mu(N;k), exactly."""
    n = np.arange(N, dtype=np.int64)
    f = np.zeros(N, dtype=np.float64)
    f[1:] = lam[1:N] * mu[(N - n)[1:]]
    C = float(f.sum())
    bad = set(prime_factors(N, spf))
    tot = 0.0
    for k in range(2, K):
        if mu[k] == 0:
            continue
        v, phi, ok = k, 1, True
        while v > 1:
            p = int(spf[v])
            if p in bad:
                ok = False
                break
            phi *= (p - 1)
            while v % p == 0:
                v //= p
        if not ok:
            continue
        r = N % k
        inner = float(f[r::k].sum()) if r else float(f[k::k].sum())
        tot += mu[k] * math.log(k) * (inner - C / phi)
    return tot, C


def main():
    global TWIN
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    TWIN = twin_constant()

    LIM = 4_000_000
    say("computing sieves to %d ..." % LIM)
    mu, _, spf = sieve(LIM)
    lam = von_mangoldt(LIM, spf)

    say()
    say("P1 / P2  --  the two constants of Section {#sec:C}")
    say("=" * 68)
    ok1 = ok2 = None
    for N in (10 ** 6, 4 * 10 ** 6):
        if N > LIM:
            continue
        pf = prime_factors(N, spf)
        S = singular_series(N, pf)
        A = A_of_N(N, pf)
        K = int(N ** THETA)
        say()
        say("N = %d   prime factors %s   S(N) = %.6f   A(N) = %.6f"
            % (N, pf, S, A))
        say("  x            A(N)*Gtilde(x)      target -S(N)   target +S(N)")
        gm = G_tilde(min(N, LIM), mu, spf, pf)
        last = None
        for x in sorted(gm):
            val = A * gm[x]
            say("  %-12d %+.6f            %+.6f      %+.6f"
                % (x, val, -S, S))
            last = val
        # The refutation rule for P1 has two clauses joined by "or", and
        # only the second was ever scored.  The first -- positive at
        # every x tested -- is restored here, and the verdict is carried
        # across the field instead of being overwritten by the last N:
        # "the largest x" is registered, "the last N" is not.
        allpos = all(A * gm[x] > 0 for x in gm)
        this1 = (not allpos) and abs(last + S) < abs(last - S)
        ok1 = this1 if ok1 is None else (ok1 and this1)
        say("  -> |A*Gt + S| = %.6f   |A*Gt - S| = %.6f   P1 %s"
            % (abs(last + S), abs(last - S),
               "HOLDS" if this1 else "REFUTED"))

        say("  K            B_log(K)            target -S(N)")
        bm = B_log(K, mu, spf, pf)
        lastb = None
        for k in sorted(bm):
            say("  %-12d %+.6f            %+.6f" % (k, bm[k], -S))
            lastb = bm[k]
        # P2 is registered as "REFUTED if positive at every K tested",
        # so the last K alone is not the registered quantity either.
        # The implementation stays stricter than the registration -- a
        # single non-negative value at the largest K refutes -- because
        # a stricter test that passes carries the looser one with it;
        # what is fixed is the field, not the threshold.
        this2 = lastb < 0 and not all(bm[k] > 0 for k in bm)
        ok2 = this2 if ok2 is None else (ok2 and this2)
        say("  -> B_log(K) is %s   P2 %s"
            % ("negative" if lastb < 0 else "positive",
               "HOLDS" if this2 else "REFUTED"))

    say()
    say("P3  --  brute-force E_3 against the two candidate identities")
    say("=" * 68)
    say("  N          E_3            rt-S(N-C)      rt+S*N+S*C     C(N)")
    ok3 = None
    for N in (200_000, 400_000, 800_000, 1_600_000):
        if N > LIM:
            continue
        pf = prime_factors(N, spf)
        S = singular_series(N, pf)
        K = int(N ** THETA)
        E3, C = E3_bruteforce(N, mu, lam, spf, K)
        n = np.arange(N, dtype=np.int64)
        rt = float((lam[1:N] * lam[(N - n)[1:]]).sum())
        cand_paper_thm = rt - S * (N - C)
        cand_paper_sec = rt + S * N + S * C
        say("  %-10d %+.4e    %+.4e    %+.4e    %+.4e"
            % (N, E3, cand_paper_thm, cand_paper_sec, C))
        ok3 = abs(E3 - cand_paper_thm) < abs(E3 - cand_paper_sec)
    say("  -> P3 %s" % ("HOLDS" if ok3 else "REFUTED"))

    say()
    say("=" * 68)
    verdict = ok1 and ok2 and ok3
    say("P1 %s   P2 %s   P3 %s"
        % ("hold" if ok1 else "REFUTED",
           "hold" if ok2 else "REFUTED",
           "hold" if ok3 else "REFUTED"))
    # The headline printed A(N)*Gtilde(1) = S(N) until 2026-08-31, without
    # the minus.  The table under it was right and the paper is right --
    # A(N)*Gtilde(x) = -1.760250 against target -S(N) = -1.760432 -- so a
    # reader who took the last line for the finding read the sign of the
    # identity backwards.  The sixth blind pass on P1 caught it.
    say("finding stands: Section {#sec:C} item (i) has A(N)*Gtilde(1) = -S(N)"
        if verdict else
        "finding refuted: A(N)*Gtilde(1) is not -S(N) on this field")

    head = [
        "STATISTIC: (a) A(N) * sum_{m<=x} mu(m)lambda(m)1_{(m,N)=1}log m / m,",
        "           (b) sum_{k<K,(k,N)=1} mu(k) log k / phi(k),",
        "           (c) E_3(alpha) by direct enumeration over k<K=N^0.56,",
        "           each against -S(N), and against the two candidate",
        "           right-hand sides of Theorem {#thm:C}.",
        "FIELD: (a),(b) N = 10^6 and 4*10^6, truncations x,K in geometric",
        "       checkpoints up to 4*10^6 and N^0.56; (c) N = 2*10^5, 4*10^5,",
        "       8*10^5, 1.6*10^6, exact integer-sieved Lambda and mu.",
        "",
    ]
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head + lines) + "\n")
    print("\nwrote %s" % OUT)
    if not verdict:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
