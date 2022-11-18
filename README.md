
⚠ WARNING: This package is very much a work in progress

## pyby

`pyby` is a re-write of the [`{BFDA}` `R`
Package](https://github.com/nicebread/BFDA).

In the future this library will be extended into a full package, but the
current version is just intended as a demonstration of how much compute
time can be reduced through the use of better algorithms and a switch in
language.

### Background

On my blog I outlined a functional re-write of the BFDA package. In
[Part I](https://blog.colling.net.nz/posts/speeding-up-bfda-part1) and
[Part II](https://blog.colling.net.nz/posts/speeding-up-bfda-part2) I
demonstrated how adopting a more functional coding style had a duel
benefit: First, it resulted in more readable code. And Second, it made
it easier to identify computational inefficiencies.

Using these principles I re-wrote the `{BFDA}` package in `R`. This led
to a large reduction in the time needed to perform the simulations (see
table below).

| VERSION                 | RUNNING TIME              |
| ----------------------- | ------------------------- |
| BFDA (original)         | 397.119s (\~6.62 minutes) |
| tidyverse-based rewrite | 95.661s (\~1.59 minutes)  |
| base R-based rewrite    | 77.866s (\~1.3 minutes)   |

This package demonstrates a re-write in Rust with an even large decrease
in compute time. For same model as outlined in the blog post
in the Python/Rust version completes the simulations
in 1.6s. The results are shown below.

```python
import pyby
import time
import pandas as pd

tic = time.perf_counter()

df = pyby.bf_sim(
    0.5,
    "t.paired",
    {"family": "Cauchy", "params": [0, 0.707], "alternative": "two.sided"},
    sampling_rule={"n_min": 50, "n_max": 300, "step_size": 5},
    alternative="two.sided",
    reps=1000,
    seed=1600,
)

df = pd.DataFrame.from_records(
    df, columns=["id", "true.ES", "n", "logBF", "emp.ES", "statistic", "p.value"]
)
toc = time.perf_counter()
print(f"Simulation performed in {toc - tic:0.4f} seconds")
print(df)
```

```md
Simulation performed in 1.6360 seconds
|       | id   | true.ES  | n    | logBF    | emp.ES    | statistic | p.value      |
|-------|------|--------- |------|-----------|----------|-----------|--------------|
| 0     |  0   | 0.5      |  50  |  5.043586 | 0.580265 | 4.103096  | 1.535430e-04 |
| 1     |  0   | 0.5      |  55  |  4.275476 | 0.515145 | 3.820419  | 3.461539e-04 |
| 2     |  0   | 0.5      |  60  |  5.397193 | 0.539404 | 4.178204  | 9.844347e-05 |
| 3     |  0   | 0.5      |  65  |  3.790283 | 0.450784 | 3.634336  | 5.565287e-04 |
| 4     |  0   | 0.5      |  70  |  3.393479 | 0.417477 | 3.492866  | 8.385545e-04 |
| ...   |...   | ...      | ...  |       ... |      ... |      ...  |          ... |
| 50995 |  999 | 0.5      | 280  | 29.076271 | 0.507505 | 8.492179  | 1.224343e-15 |
| 50996 |  999 | 0.5      | 285  | 28.385162 | 0.496423 | 8.380589  | 2.467588e-15 |
| 50997 |  999 | 0.5      | 290  | 27.564108 | 0.484431 | 8.249554  | 5.680903e-15 |
| 50998 |  999 | 0.5      | 295  | 26.432018 | 0.469911 | 8.070985  | 1.798275e-14 |
| 50999 |  999 | 0.5      | 300  | 26.020612 | 0.461992 | 8.001939  | 2.724417e-14 |

[51000 rows x 7 columns]
```
