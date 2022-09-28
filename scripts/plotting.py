import json

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def readpositions(filename, axes):
    with open(filename, 'r') as f:
        data = json.load(f)
    positions = list(map(lambda x: x["position"], data))
    if axes == "all":
        return positions
    elif axes == "x":
        return list(map(lambda x: x[0], positions))
    elif axes == "y":
        return list(map(lambda x: x[1], positions))
    elif axes == "z":
        return list(map(lambda x: x[2], positions))
    else:
        pass


def plot(files):
    fig = plt.figure(figsize=(16, 10))
    ax = plt.axes(projection="3d")
    for file in files:
        x = readpositions(file, "x")
        y = readpositions(file, "y")
        z = readpositions(file, "z")
        ax.scatter3D(x, y, z, label=file)
    ax.legend(loc="best")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    # fig.savefig("images.pdf")
    return ax


if __name__ == "__main__":
    a = plot(["../x.json", "../y.json"])
