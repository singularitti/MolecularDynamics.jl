import json

import matplotlib.pyplot as plt

params = {
    "figure.figsize": (8, 5),
    "figure.autolayout": True,
    "figure.dpi": 200,
    "legend.fontsize": "small",
    "text.usetex": True,
}
plt.rcParams.update(params)


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


def plot(files, ax=plt.axes(projection="3d")):
    fig = plt.figure(figsize=(16, 10))
    for file in files:
        x = readpositions(file, "x")
        y = readpositions(file, "y")
        z = readpositions(file, "z")
        ax.scatter3D(x, y, z, label=file)
    # See https://stackoverflow.com/a/63625222/3260253
    ax.set_box_aspect([1, 1, 1])
    ax.legend(loc="best")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    fig.savefig("images.pdf")
    return fig, ax


if __name__ == "__main__":
    data = readpositions("../p.json", "all")
    arr = ["../p.json"]
    # r = range(10, 21, 10)
    r = (29, 49)
    for i in r:
        arr.append("../" + str(i+1) + ".json")
    fig, ax = plot(arr)
    for i in r:
        ax.scatter(*data[i], label=str(i+1), s=40)
    ax.legend(loc="best")
