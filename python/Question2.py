import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from libraryWrapper import WorkerLibraryWrapper
from Analytical import cylindricalInclusion


def main():
    inclusionModuli = np.array([1e-15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    numerical = np.array([])
    analytical = np.array([])
    wrapper = WorkerLibraryWrapper()
    wrapper.verbose = False
    for E in tqdm(inclusionModuli):
        wrapper.inclusionModulus = E
        result = wrapper.worker()
        numerical = np.append(numerical, result)
        analytical = np.append(
            analytical,
            cylindricalInclusion(
                E_m=wrapper.matrixModulus,
                nu_m=wrapper.matrixPoisson,
                E_i=wrapper.inclusionModulus,
                nu_i=wrapper.inclusionPoisson,
            ),
        )
    # Plot
    plt.rcParams["font.size"] = 24
    plt.figure(figsize=(10, 8))
    plt.gca().set_aspect("equal", adjustable="box")
    params = {
        "markersize": 12,
        "markeredgewidth": 2,
        "markerfacecolor": "none",
    }
    plt.plot(
        inclusionModuli,
        analytical,
        label="Analytical",
        marker="o",
        **params,
    )
    plt.plot(
        inclusionModuli,
        numerical,
        label="Numerical",
        marker="x",
        **params,
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
