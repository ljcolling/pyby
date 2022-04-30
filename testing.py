from pyby import gen_data as r
import pandas as pd

output = r(20, 10)


# def give_index(index, group):
#     if group == 2:
#         return index - 10
#     return index


# output = [
#     (give_index(index + 1, group), group, value)
#     for index, (group, value) in enumerate(output)
# ]

pd.DataFrame.from_records(output, columns=["n", "val", "running_mean"]).to_csv(
    "data2.csv", index=False
)
