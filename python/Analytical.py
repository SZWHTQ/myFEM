import numpy as np
import matplotlib.pyplot as plt


def getLame(E, nu):
    mu = E / (2 * (1 + nu))
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    return lam, mu


def cylindricalInclusion(E_m=1, nu_m=0.3, E_i=1, nu_i=0.2, sigma_0=1, a=1):
    _, mu_m = getLame(E_m, nu_m)
    _, mu_i = getLame(E_i, nu_i)
    u_r = (
        sigma_0 * a * (1 - nu_m) / mu_m * (1 - 2 * nu_i) / (1 - 2 * nu_i + mu_i / mu_m)
    )
    deltaU = (
        np.pi
        * sigma_0
        * a
        * u_r
        * (
            1
            - (1 + nu_m) * E_i / (1 + nu_i) / E_m
            - 2 * (nu_i - nu_m) * (1 + nu_m) * E_i / (1 + nu_i) / (1 - 2 * nu_i) / E_m
        )
    )
    return -deltaU


def circle(E_m=1, nu_m=0.3, E_i=1, nu_i=0.2, sigma_0=1, a=1):
    L_m, mu_m = getLame(E_m, nu_m)
    L_i, mu_i = getLame(E_i, nu_i)

    A = 0.5 * sigma_0 / (L_m + mu_m)
    B = (
        0.5
        * sigma_0
        * a
        * a
        * (L_m + mu_m - L_i - mu_i)
        / ((L_m + mu_m) * (L_i + mu_i + mu_m))
    )

    deltaU = -1 * np.pi * a * sigma_0 * (A * a + B / a) * (E_m - E_i) / E_m
    return -deltaU


def main():
    E = [1e-15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    Y1 = [cylindricalInclusion(E_i=E_i) for E_i in E]
    Y2 = [
        -5.930960,
        -0.568125,
        0.549193,
        1.032815,
        1.302745,
        1.474987,
        1.594445,
        1.682160,
        1.749299,
        1.802343,
        1.845309,
    ]
    ## Plot
    # Settings
    plt.rcParams["font.size"] = 24
    plt.figure(figsize=(10, 8))
    plt.gca().set_aspect("equal", adjustable="box")
    # Curves
    plt.plot(
        E,
        Y1,
        label="Analytical",
        marker="o",
        markersize=12,
        markeredgewidth=2,
        markerfacecolor="none",
        linestyle="-",
        linewidth=2,
        color="blue",
    )
    plt.plot(
        E,
        Y2,
        label="Numerical",
        marker="x",
        markersize=12,
        markeredgewidth=2,
        linestyle="-",
        linewidth=2,
        color="red",
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
