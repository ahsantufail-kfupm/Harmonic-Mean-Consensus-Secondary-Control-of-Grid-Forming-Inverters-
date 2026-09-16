"""Verification script for Remark (Weight magnitude versus weight shape).

Checks, over 500 random connected graphs, whether the harmonic-mean weighting
retains the largest lambda2 when every power-mean weighting is renormalised to
the same total edge-weight budget sum_e a_e = |E|.

Reference run (numpy default_rng(seed=7), n in [6,25), G(n,p) with p in [0.15,0.5)):
  RAW (as computed, no renormalisation): harmonic largest 500/500, full ordering 500/500
  NORMALISED (equal budget):             harmonic largest 384/500
                                         harmonic/family-best ratio mean 0.9906, min 0.7960
Re-run this (or port to MATLAB) and confirm before the numbers in the Remark are kept.
"""
import numpy as np, networkx as nx

rng = np.random.default_rng(7)
def lam2(L): return np.linalg.eigvalsh(L)[1]
def wlap(G, w):
    n = G.number_of_nodes(); L = np.zeros((n, n))
    for (i, j), we in zip(G.edges(), w):
        L[i,i]+=we; L[j,j]+=we; L[i,j]-=we; L[j,i]-=we
    return L
def pm(di, dj, p):
    return np.sqrt(di*dj) if p == 0 else (0.5*(di**p + dj**p))**(1.0/p)

ps = [-1, 0, 1, 2]; N = 500
raw_wins = norm_wins = raw_ord = norm_ord = 0; ratios = []
for _ in range(N):
    n = rng.integers(6, 25)
    while True:
        G = nx.gnp_random_graph(int(n), rng.uniform(0.15, 0.5), seed=int(rng.integers(1e9)))
        if nx.is_connected(G): break
    G = nx.convert_node_labels_to_integers(G)
    deg = dict(G.degree()); m = G.number_of_edges()
    lr, ln = [], []
    for p in ps:
        w = np.array([1.0/pm(deg[i], deg[j], p) for i, j in G.edges()])
        lr.append(lam2(wlap(G, w)))
        ln.append(lam2(wlap(G, w*(m/w.sum()))))
    raw_wins  += lr[0] == max(lr);  norm_wins += ln[0] == max(ln)
    raw_ord   += all(lr[k] >= lr[k+1]-1e-12 for k in range(3))
    norm_ord  += all(ln[k] >= ln[k+1]-1e-12 for k in range(3))
    ratios.append(ln[0]/max(ln))
print(f"RAW : harmonic largest {raw_wins}/{N}, ordering {raw_ord}/{N}")
print(f"NORM: harmonic largest {norm_wins}/{N}, ordering {norm_ord}/{N}")
print(f"NORM harmonic/best ratio: mean {np.mean(ratios):.4f}, min {np.min(ratios):.4f}")
