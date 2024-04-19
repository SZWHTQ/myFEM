"""Question3_2.py, a/b=1/3"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from tqdm import tqdm
from libraryWrapper import WorkerLibraryWrapper
import scipy as sp


def getCurve(Ksi, inclusionModulus):
    result = np.array([])
    wrapper = WorkerLibraryWrapper()
    wrapper.inclusionModulus = inclusionModulus
    wrapper.a_b = 1 / 3
    wrapper.meshSize = 0.5
    wrapper.refinementFactor = 20
    # wrapper.verbose = True
    for xi in Ksi:
        wrapper.ksi = xi
        result = np.append(result, wrapper.run())
    return result


def main():
    plt.rcParams["font.size"] = 24
    plt.figure(figsize=(16, 9))
    # plt.set_aspect("equal", adjustable="box")

    params = {
        "linewidth": 3,
        "markersize": 16,
        "markeredgewidth": 3,
        "markerfacecolor": "none",
    }
    colors = ["red", "green", "blue", "black"]
    lineStyles = ["--", "-.", ":"]

    Ksi = np.arange(7, 24, 1, dtype=np.float64)
    inclusionModulus = np.arange(0, 11, 1, dtype=np.float64)
    inclusionModulus[0] = 1e-15
    Curves = []
    for E in tqdm(inclusionModulus):
        Curves.append(getCurve(Ksi, E))

    sp.io.savemat(
        "Question3_2.mat",
        {"Ksi": Ksi, "Curves": Curves, "inclusionModulus": inclusionModulus},
    )

    ax = plt.gca()
    for i, curve in enumerate(Curves):
        ax.plot(
            Ksi,
            curve,
            label=f"{inclusionModulus[i]}",
            **params,
            color=colors[i % len(colors)],
            linestyle=lineStyles[i % len(lineStyles)],
        )

    # Set the first y-axis label and scale
    ax.set_xlabel(r"$\xi$")
    ax.set_ylabel(r"$\Delta U$")
    xLimit = [6.5, 23.5]
    yLimit = [-165, 20]
    xRange = [7, 24]
    yRange = [-160, 21]
    ax.set_xlim(*xLimit)
    ax.set_ylim(*yLimit)
    ax.set_xticks(np.arange(*xRange, 2))
    ax.set_yticks(np.arange(*yRange, 30))
    ax.minorticks_on()
    ax.set_xticks(np.arange(*xRange, 1), minor=True)
    ax.set_yticks(np.arange(*yRange, 15), minor=True)

    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.legend(loc="lower right")
    # Set legend font size
    for text in ax.get_legend().get_texts():
        text.set_fontsize(16)

    # Add inset of the same curve in the lower left corner
    axins = inset_axes(
        ax,
        width="50%",
        height="50%",
        loc="lower left",
        bbox_to_anchor=(0.25, 0.05, 1, 1),
        bbox_transform=ax.transAxes,
    )

    # Plot the same curve on the overlapping subgraph
    for i, curve in enumerate(Curves):
        axins.plot(
            Ksi,
            curve,
            **params,
            color=colors[i % len(colors)],
            linestyle=lineStyles[i % len(lineStyles)],
        )

    xLimitIns = [6.9, 9.1]
    yLimitIns = [1.9, 12.1]
    xRangeIns = [7, 9.1]
    yRangeIns = [2, 12]

    # 绘制方框，标注出原图中对应子图的区域
    # 这里的(2, -10)是子图的左下角坐标，(6-2, 20-(-10))是宽度和高度
    rect = Rectangle(
        (xLimitIns[0], yLimitIns[0]),
        xLimitIns[1] - xLimitIns[0],
        yLimitIns[1] - yLimitIns[0],
        linewidth=2,
        edgecolor="r",
        facecolor="none",
    )
    ax.add_patch(rect)

    # 在原图中添加箭头，指向重叠子图
    # xy参数是箭头指向的位置，xytext是文本的位置，这里我们没有添加文本，只有箭头
    ax.annotate(
        "",
        xy=(
            xLimit[0] + (xLimit[1] - xLimit[0]) * 0.27,
            yLimit[0] + (yLimit[1] - yLimit[0]) * 0.57,
        ),
        xytext=(xRangeIns[1], yRangeIns[0]),
        arrowprops=dict(
            color="red", arrowstyle="->,head_width=0.4,head_length=1", linewidth=2
        ),
    )

    axins.set_xlim(*xLimitIns)
    axins.set_ylim(*yLimitIns)
    axins.set_xticks(np.arange(*xRangeIns, 0.25))
    axins.set_yticks(np.arange(*yRangeIns, 1))
    axins.minorticks_on()
    axins.set_xticks(np.arange(*xRangeIns, 0.5), minor=True)
    axins.set_yticks(np.arange(*yRangeIns, 2), minor=True)
    axins.grid(True, which="both", linestyle="--", linewidth=0.25)
    # Make it look clearer
    axins.tick_params(axis="both", which="major", labelsize=8)
    axins.tick_params(axis="both", which="minor", labelsize=6)

    plt.show()


if __name__ == "__main__":
    main()
