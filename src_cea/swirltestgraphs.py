import numpy as np
import matplotlib.pyplot as plt

# Define the function
def y(x, A):
    return 1 / np.sqrt(A**2/(1-x) + 1/x**2)

# Define the range of x values to plot
x = np.linspace(0.01, 0.99, 1000)

# Define the values of A to plot
A_values = [2.0, 0.8, 0.4, 0.2]

# Create a figure and axes object
fig, ax = plt.subplots()

# Plot the function for each value of A
for A in A_values:
    ax.plot(x, y(x, A), label=f"A={A}")
ax.plot(x, x*np.sqrt(x/(2-x)), label=f"ideal")

# Set the title and axis labels
ax.set_title("Graph of y = 1/(A^2/(1-x)+1/x^2)")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Add a legend and gridlines
ax.legend()
ax.grid()

# Show the plot
#plt.show()

# Define the function
def phi(x, y):
    return y * (y / (2 - y)) ** 0.5 - 1 / ((x**2 / (1-y) + 1/y**2) ** 0.5)

# Define the ranges of x and y values to plot
x = np.linspace(0.01, 11.99, 1000)
y = np.linspace(0.01, 0.99, 1000)

# Create a meshgrid of x and y values
X, Y = np.meshgrid(x, y)

# Compute the function values for the meshgrid
Z = phi(X, Y)
# Create a figure and axes object
fig, ax = plt.subplots()

# Plot the contour lines where f(x, y) = 0
ax.contour(X, Y, Z, levels=[0], colors='black')

# Set the title and axis labels
ax.set_title("Graph of y*(y/(2-y))^0.5 = 1/(x^2/(1-y) + 1/y^2)^0.5")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Show the plot
#plt.show()
#--------------------------------------------------------------------------------------------------------------------------------------
# Define the functions
def f(φ, A):
    return φ * (φ / (2 - φ)) ** 0.5 - 1 / ((A ** 2 / (1 - φ) + 1 / φ ** 2) ** 0.5)

def μ(φ):
    return φ * (φ / (2 - φ)) ** 0.5

def U(φ):
    return (2 * (1 - φ) * (2 - φ)) ** 0.5

def V(φ):
    return (φ / (2 - φ)) ** 0.5

# Define the range of A values to plot
A = np.linspace(0.01, 2, 1000)

# Define the range of φ values to plot
φ = np.linspace(0.01, 1.99, 1000)

# Create a figure and axes object
fig, axs = plt.subplots(2, 2)

# Plot the first function on the first axis
axs[0, 0].plot(A, f(φ=0.5, A=A))
axs[0, 0].set_title("φ = 0.5")
axs[0, 0].set_xlabel("A")
axs[0, 0].set_ylabel("y")

# Plot the second function on the second axis
axs[0, 1].plot(φ, μ(φ=φ))
axs[0, 1].set_title("μ(φ)")
axs[0, 1].set_xlabel("φ")
axs[0, 1].set_ylabel("μ")

# Plot the third function on the third axis
axs[1, 0].plot(φ, U(φ=φ))
axs[1, 0].set_title("U(φ)")
axs[1, 0].set_xlabel("φ")
axs[1, 0].set_ylabel("U")

# Plot the fourth function on the fourth axis
axs[1, 1].plot(φ, V(φ=φ))
axs[1, 1].set_title("V(φ)")
axs[1, 1].set_xlabel("φ")
axs[1, 1].set_ylabel("V")

# Adjust the layout and spacing of the subplots
fig.tight_layout(pad=1.0)

# Show the plot
plt.show()