import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from libraryWrapper import WorkerLibraryWrapper
from Analytical import cylindricalInclusion


# inclusionModuli = np.array([1e-15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
inclusionModuli = np.arange(0, 11, 1, dtype=np.float64)
inclusionModuli[0] = 1e-15
numerical = np.array([])
analytical = np.array([])
wrapper = WorkerLibraryWrapper()
wrapper.meshSize = 2
wrapper.refinementFactor = 10
for E in tqdm(inclusionModuli):
    wrapper.inclusionModulus = E
    numerical = np.append(numerical, wrapper.run())
    analytical = np.append(
        analytical,
        cylindricalInclusion(
            E_m=wrapper.matrixModulus,
            nu_m=wrapper.matrixPoisson,
            E_i=wrapper.inclusionModulus,
            nu_i=wrapper.inclusionPoisson,
        ),
    )
relativeError = np.abs(analytical - numerical) / np.abs(analytical) * 100

# Plot
plt.rcParams["font.size"] = 24
plt.figure(figsize=(16, 9))
ax1 = plt.gca()
# ax1.set_aspect("equal", adjustable="box")

params = {
    "linewidth": 3,
    "markersize": 16,
    "markeredgewidth": 3,
    "markerfacecolor": "none",
}

ax1.plot(
    inclusionModuli,
    analytical,
    label="Analytical",
    marker="o",
    **params,
)

ax1.plot(
    inclusionModuli,
    numerical,
    label="Numerical",
    marker="x",
    **params,
)

# Set the first y-axis label and scale
ax1.set_xlabel(r"$E_I/E_M$")
ax1.set_ylabel(r"$\Delta U$")
ax1.set_ylim(-6.5, 2.5)
ax1.set_xticks(np.arange(0, 11, 2))
ax1.set_yticks(np.arange(-6, 3, 2))
ax1.minorticks_on()
ax1.set_xticks(np.arange(0, 11, 1), minor=True)
ax1.set_yticks(np.arange(-6, 3, 1), minor=True)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5)
ax1.set_axisbelow(True)

# Create the second y-axis
ax2 = ax1.twinx()
ax2.plot(
    inclusionModuli,
    relativeError,
    label="Relative Error",
    color="red",
    marker="s",
    **params,
)

# Set the second y-axis label and scale
ax2.set_ylabel("Relative Error(%)")
ax2.set_yticks(np.linspace(0, 4, 5))
ax2.set_ylim(-0.25, 4.25)

# Legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc="lower right")

plt.show()
