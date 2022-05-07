import pyby
import time
import pandas as pd

tic = time.perf_counter()
df = pyby.bf_sim(
    0.5,
    "paired_t",
    {"name": "cauchy", "params": [0, 0.707]},
    10,
    300,
    10,
    "two.sided",
    1000,
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
