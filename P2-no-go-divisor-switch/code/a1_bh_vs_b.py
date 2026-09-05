# -*- coding: utf-8 -*-
"""
a1_bh_vs_b.py  --  pass4, blind mathematical re-verification.

PRE-REGISTRATION (fixed before the run).

WHAT IS MEASURED.  At N in {2e5,4e5,8e5,1.6e6,3.2e6}, theta'=0.56,
K=floor(N**0.56), over moduli k<K with (k,N)=1 (and separately over the
squarefree such k):

  A(N;k)    = sum_{n<N, n = N (mod k)} Lambda(n) mu(N-n)
  C(N)      = sum_{n<N} Lambda(n) mu(N-n)
  E_mu(N;k) = A(N;k) - C(N)/phi(k)
  H(N;k)    = sum_{1<=m<N/k, (m,k)=1} Lambda(N-mk) mu(m)
  B_all  = sum_k (log k)|E_mu|,  B_sqf  = same over squarefree k
  BH_all = sum_k (log k)|H|,     BH_sqf = same over squarefree k

and the two ratios P2's Measurement 11 prints, with
S(N) the singular series and A(N)=prod_{q not| N}(1-1/(q(q-1))).

WHAT WOULD FALSIFY THE FINDING UNDER TEST.  The finding is: P2
Measurement 11 writes BOTH of its ratio lists with the numerator
B_H(N)=sum(log k)|H(N;k)|, but its second ("old") list is in fact
built from B(N)=sum(log k)|E_mu(N;k)| -- P1's eq:(9) quantity, a
different sum.  FALSIFIED if BH/(S(1-A)N) reproduces the printed
2.1591, 1.9747, 1.9500, 1.7483, 1.5798 to the digits printed.
CONFIRMED if it does not and B/(S(1-A)N) does.

PREDICTION.  E_mu(N;k) = mu(k)H(N;k) - C(N)/phi(k); inside an absolute
value the subtracted mean term does not cancel, so B != B_H by a
quantity of order |C(N)|(log K)^2, i.e. of order 1e-2 N or more here,
which is larger than the gap between the two printed lists.

NULL.  None applies and none would mean anything: every quantity here
is a deterministic arithmetic sum with no sign input of its own, and a
sign control would alter both sides of every comparison alike.  The
reference is the printed digit string, which is what the test compares
against.
"""
import numpy as np, math, sys

def sieve_all(n):
    sieve = np.ones(n+1, dtype=bool); sieve[:2] = False
    for i in range(2, int(n**0.5)+1):
        if sieve[i]: sieve[i*i::i] = False
    primes = np.nonzero(sieve)[0]
    mu = np.ones(n+1, dtype=np.int64); mu[0] = 0
    phi = np.arange(n+1, dtype=np.int64)
    lam = np.zeros(n+1, dtype=np.float64)
    for p in primes:
        p = int(p)
        mu[p::p] *= -1
        phi[p::p] -= phi[p::p]//p
        if p*p <= n: mu[p*p::p*p] = 0
        q = p
        lp = math.log(p)
        while q <= n:
            lam[q] = lp; q *= p
    return primes, mu, phi, lam

def factor(n, primes):
    f = []; m = n
    for p in primes:
        p = int(p)
        if p*p > m: break
        if m % p == 0:
            f.append(p)
            while m % p == 0: m //= p
    if m > 1: f.append(m)
    return f

def main():
    Ns = [200000, 400000, 800000, 1600000, 3200000]
    LIM = max(Ns)
    sys.stderr.write("sieving to %d ...\n" % LIM)
    primes, mu, phi, lam = sieve_all(LIM)
    print(__doc__.strip())
    print()
    print("  %-9s %-6s %-9s %-9s %-9s %-9s %-9s %-9s" %
          ("N","K","B_all/N","B_sqf/N","BH_all/N","BH_sqf/N","S(N)","1-A(N)"))
    print("  "+"-"*82)
    rows = []
    for N in Ns:
        K = int(N**0.56)
        n = np.arange(1, N)
        Cn = float(np.dot(lam[1:N], mu[N-n]))
        pf = factor(N, primes)
        Ball=Bsqf=BHall=BHsqf=0.0
        for k in range(2, K):
            if math.gcd(k, N) != 1: continue
            M = (N-1)//k
            ms = np.arange(1, M+1)
            u = ms*k
            nn = N - u
            Ak = float(np.dot(lam[nn], mu[u]))
            keep = np.ones(M, dtype=bool)
            for p in factor(k, primes):
                keep[p-1::p] = False
            Hk = float(np.dot(lam[nn][keep], mu[ms][keep]))
            Ek = Ak - Cn/phi[k]
            lg = math.log(k)
            Ball += lg*abs(Ek); BHall += lg*abs(Hk)
            if mu[k] != 0:
                Bsqf += lg*abs(Ek); BHsqf += lg*abs(Hk)
        S = 2.0; A = 1.0
        for q in primes:
            q = int(q)
            if q > 2: S *= (1.0 - 1.0/float((q-1)**2))
            if N % q != 0: A *= (1.0 - 1.0/float(q*(q-1)))
        for p in pf:
            if p > 2: S *= (1.0 + 1.0/float(p-2))
        print("  %-9d %-6d %-9.4f %-9.4f %-9.4f %-9.4f %-9.6f %-9.6f" %
              (N,K,Ball/N,Bsqf/N,BHall/N,BHsqf/N,S,1-A))
        rows.append((N,K,Ball,Bsqf,BHall,BHsqf,S,A,Cn))
        sys.stderr.write("  done N=%d\n" % N)
    print()
    print("  P2 Measurement 11's two lists, recomputed four ways")
    print("  %-9s %-11s %-8s | %-13s %-13s %-13s %-8s" %
          ("N","BH_sqf/(SN)","pub new","B_all/(S(1-A)N)","B_sqf/(...)","BH_sqf/(...)","pub old"))
    pubnew = [0.4578,0.4064,0.4079,0.3769,0.3338]
    pubold = [2.1591,1.9747,1.9500,1.7483,1.5798]
    for i,(N,K,Ball,Bsqf,BHall,BHsqf,S,A,Cn) in enumerate(rows):
        d1 = S*N; d2 = S*(1-A)*N
        print("  %-9d %-11.4f %-8.4f | %-13.4f %-13.4f %-13.4f %-8.4f" %
              (N, BHsqf/d1, pubnew[i], Ball/d2, Bsqf/d2, BHsqf/d2, pubold[i]))
    print()
    print("  the size of the gap the paper's formula hides")
    for (N,K,Ball,Bsqf,BHall,BHsqf,S,A,Cn) in rows:
        print("  N=%-9d |C|/N=%.6f  (B_all-BH_sqf)/N=%+.4f  (B_sqf-BH_sqf)/N=%+.4f"
              % (N, abs(Cn)/N, (Ball-BHsqf)/N, (Bsqf-BHsqf)/N))

main()
