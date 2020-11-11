import numpy as np
import pandas as pd


def de_dup(lists):
    _set =  [tuple(x) for x in lists]
    print(type(_set))
    _set = set(_set)
    res = np.array(list(_set))
    return res


if __name__ == "__main__":
    # data = np.loadtxt("pareto.txt")
    print("load start!")
    # data = np.loadtxt("const_opt_result.csv",delimiter=",")[:, 1:]
    data = pd.read_csv("const_opt_result.csv", skiprows=1, header=None).values[:, 1:]
    print("de_duplication start!")
    data = de_dup(data)
    data = data[np.argsort(data[:, 0])]

    # np.savetxt("pareto_modified.txt", data, header=" ", footer=" ", fmt="%12.7f")
    np.savetxt("duped.csv", data, delimiter=",", header=" ", footer=" ", fmt="%.7f")


