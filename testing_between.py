from pyby import bf_sim_independent
from math import sqrt

def eff_size(m1, m2, s1, s2, n):
    n1 = n
    n2 = n
    md_diff = m1 - m2
    sd_pooled = sqrt((((n1 - 1) * s1**2.0) + ((n2 - 1) * s2**2.0)) / (n1 + n2 - 2))
    d = md_diff / sd_pooled
    return d




sum(
    [
        eff_size(m1, m2, s1, s2, n)
        for g, n, v1, v2, m1, m2, s1, s2 in bf_sim_independent(0.5, 50, 50, 1, 1000)
    ]
) / 1000
