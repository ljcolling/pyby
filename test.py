import pyby
import time
import pandas as pd

# print(a)
tic = time.perf_counter()
max_n = 200
batches = 1000

# max_n = 20
# batches = 10
effsize = 0.5
a = pyby.gen_data(effsize, batches * max_n, max_n, 69)
# a = pyby.gen_data(200, 1000, 69)
# 200 * 100000

toc = time.perf_counter()
print(f"Simulation performed in {toc - tic:0.4f} seconds")
a = pd.DataFrame.from_records(
    a, columns=["id", "true.ES", "n", "logBF", "emp.ES", "statistic", "p.value"]
)
toc = time.perf_counter()
print(f"Simulation performed in {toc - tic:0.4f} seconds")
print(a)
a.to_csv("./testdata2.csv", index= False)
