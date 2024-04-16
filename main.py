import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from WorkerLibraryWrapper import WorkerLibraryWrapper
from Circle import theoretical


def main():
    Es = np.array([1e-15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    W = WorkerLibraryWrapper()
    numerical = np.array([])
    analytical = np.array([theoretical(E_i=E_i) for E_i in Es])
    for E in tqdm(Es):
        W.inclusiveModulus = E
        result = W.worker()
        numerical = np.append(numerical, result)
    # Plot
    plt.rcParams["font.size"] = 24
    plt.figure(figsize=(10, 8))
    plt.gca().set_aspect("equal", adjustable="box")
    plt.plot(
        Es,
        analytical,
        label="Analytical",
        marker="o",
        markersize=12,
        markeredgewidth=2,
        markerfacecolor="none",
    )
    plt.plot(
        Es,
        numerical,
        label="Numerical",
        marker="x",
        markersize=12,
        markeredgewidth=2,
        markerfacecolor="none",
    )
    # Labels
    plt.xlabel(r"$E_I/E_M$")
    plt.ylabel(r"$\Delta U$")
    # Ticks
    plt.xticks(np.arange(0, 11, 2))
    plt.yticks(np.arange(-6, 3, 2))
    # Turn on minor ticks
    plt.minorticks_on()
    plt.xticks(np.arange(0, 11, 1), minor=True)
    plt.yticks(np.arange(-6, 3, 1), minor=True)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.gca().set_axisbelow(True)
    # Adjust grid line color
    plt.legend(loc="lower right")
    plt.show()


if __name__ == "__main__":
    main()
