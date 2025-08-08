import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

def upperbound(p2, n2, p3, n3, N, l, E, const=1): 
  """Calculates the upper bound."""
  r = min(p2**(n2-1)*(p2-1), p3**(n3-1)*(p3-1))
  return N * l * E * (4*p2*p3*r**2 + 6*p2*r*p3**(n3) + 2*r*p2 + 3*p2**n2) * const

def prob_upperbound(p2, n2, p3, n3, N, l, E, const=1):
  """Calculates the probabilistic upper bound."""
  r = min(p2**(n2-1)*(p2-1), p3**(n3-1)*(p3-1))
  return np.sqrt(N * l) * E * r**3 * const # * (4*p2*p3*r**2 + 6*p2*r*p3**(n3) + 2*r*p2 + 3*p2**n2) 

# --- Parameters ---
p1, n1 = 31, 1
p2, n2 = 3, 2
p3, n3 = 7, 1
N_calc = p1**(n1-1) * (p1-1) * p2**(n2-1) * (p2-1) * p3**(n3-1) * (p3-1)
l = 100
E = 7

# --- Line Data ---
x_line = np.array([i for i in range(1,21, 2)])
y_line = x_line * prob_upperbound(p2, n2, p3, n3, N_calc, l, E, 0.8) + E

# --- Scatter Data ---
#x_scatter = np.arange(1, 11)
#y_scatter = np.array([
#    4149.25, 5577.30, 6180.20, 7822.80, 8972.15,
#    9644.25, 11190.15, 11499.80, 13705.00, 12766.80
#])

x_scatter = np.array([i for i in range(1,21, 2)])
# y_scatter = np.array([[21124, 33018, 59926, 66043, 84767, 180318, 180367, 150189, 181938,139218]])
y_scatter = np.array([129556, 109067, 132342, 126396, 170701, 208909, 233064, 235783, 355840,220509])
# Create a figure with higher resolution
plt.figure(figsize=(10, 6), dpi=150)

# --- Plotting with Improvements ---
# 1. Plotted points and line with distinct colors and descriptive labels
plt.scatter(x_scatter, y_scatter, color="steelblue", label="Máximo da norma infinita", zorder=5)
plt.plot(x_line, y_line, color="darkorange", linestyle='--', label="Norma infinita Probabilística")

# --- Grid and Ticks ---
plt.grid(which="both", axis="y", linestyle="--", linewidth=0.5)
plt.xticks(np.arange(1, 21, 1)) # Set x-ticks to match data range 1-10
plt.grid(which="major", axis="x", linestyle="--", linewidth=0.5)

# --- Titles and Labels ---
plt.title("Norma infinita de Ruído por Profundidade")
plt.ylabel("$||e||_{\infty}$") # Corrected Y-axis label to reflect data
plt.xlabel("Profundidade")

# 2. Added a standard legend for plot elements inside the plot area
plt.legend(loc='upper left')

# --- Parameter Information Box ---
# Using a text box is a cleaner way to display static parameter information
p1, n1 = 31, 1
p2, n2 = 3, 2
p3, n3 = 7, 1
prime_powers = [(p1, n1), (p2, n2), (p3, n3)]
m1 , m2, m3 = p1**n1, p2**n2, p3**n3
n =p1**(n1-1) * (p1-1) * p2**(n2-1) * (p2-1) * p3**(n3-1) * (p3-1)
r = 6
Q = 2**100
sigma = 3.2

param_text = (
    f"Parametros\n"
    f"$q$ = {p1**n1}\n"
    f"$p_2$ = {p2}, $n_2$ = {n2}\n"
    f"$p_3$ = {p3}, $n_3$ = {n3}\n"
    f"r = {r}\n"
    f"N = {n}\n"
    f"Q = $2^{{100}}$\n"
    f"$\\sigma$ = {sigma}"
)

# Place the text box outside the plot on the right side
plt.text(1.02, 0.98, param_text, transform=plt.gca().transAxes, fontsize=9,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

# Adjust layout to prevent the text box from being cut off
plt.tight_layout(rect=[0, 0, 0.85, 1])

# Show the plot
#plt.show()
plt.savefig("bbbbb.jpg")