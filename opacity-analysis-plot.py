#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
names: "Opacity", "Surface-Temperature_Grey[K]", "Surface-Temperature_Window[K]"
"""
table = pd.read_csv("opacity-analysis.csv", delim_whitespace=True)

plt.figure()
plt.plot(table["Opacity"], table["Surface-Temperature_Grey[K]"], label="Gray", marker="1")
plt.plot(table["Opacity"], table["Surface-Temperature_Window[K]"], label="Window", marker="1")
plt.title("Equilibrium Surface Temperature for Various Opacities")
plt.xlabel("Opacity")
plt.ylabel("Surface Temperature [K]")
plt.legend()
plt.savefig("opacity analysis.pdf")
plt.show(block=False)

plt.show()
