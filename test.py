import pyby
import time
import pandas as pd

tic = time.perf_counter()

df = pyby.bf_sim(
    0.5,
    "t.paired",
    {"family": "Cauchy", "params": [0, 0.707], "alternative": "two.sided"},
    sampling_rule = {"n_min": 50, "n_max": 300, "step_size": 5},
    alternative = "two.sided",
    reps = 1000,
    seed = 1600,
)

if df is not None:
    df = pd.DataFrame.from_records(
        df, columns=["id", "true.ES", "n", "logBF", "emp.ES", "statistic", "p.value"]
    )
    toc = time.perf_counter()
    print(f"Simulation performed in {toc - tic:0.4f} seconds")
    print(df)
else:
    Exception("Error")




