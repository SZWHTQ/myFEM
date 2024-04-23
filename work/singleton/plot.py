import numpy as np
import matplotlib.pyplot as plt

# Data for plotting
x = np.array([2.1, 3, 4, 5, 6, 7, 8, 9, 10])

# # a=b=1, 1 0.3 2 0.2
# y = np.array(
#     [
#         0.129304,
#         0.133425,
#         0.135283,
#         0.136164,
#         0.136650,
#         0.136943,
#         0.137124,
#         0.137234,
#         0.137298,
#     ]
# )

# a=1 b=1/3, 1 0.3 2 0.2
y = np.array(
    [
        0.045672,
        0.046060,
        0.046240,
        0.046328,
        0.046376,
        0.046403,
        0.046418,
        0.046425,
        0.046429,
    ]
)


# plt.rcParams['text.usetex'] = True
plt.xlabel("B")
plt.ylabel(r"$\Delta U$")
plt.plot(x, y, "o-")

plt.show()
