import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

# Parameters of equation
a = 1.5
l = 2 * np.pi
T = 2
# Steps
N = 50
M = 10 + 2 * pow(N, 2)
x = np.linspace(0, l, N)
t = np.linspace(0, T, M)
h = x[1]
tau = t[1]
# Condition of stability
print(tau <= pow(h, 2) / (2 * pow(a, 2)))
# Temperature
u = []
# Initial condition
u.append([3 * (l - i) for i in x])


# Function
def f(x, t):
    return t * np.sin(x)


'''Explicit method'''
for j in range(M - 1):
    x_array = [0] * N
    for i in range(1, N - 1):
        x_array[i] = u[j][i] + (tau * a ** 2) / pow(h, 2) * (u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) + tau * f(x[i],
                                                                                                                t[j])
    x_array[0] = x_array[1]
    x_array[N - 1] = x_array[N - 2]
    u.append(x_array)

fig, ax = plt.subplots()

line, = ax.plot(x, u[0])


def animate(i):
    line.set_ydata(u[i])  # update the data
    return line,


ani = animation.FuncAnimation(fig, animate, interval=1, frames=M, repeat=False)

plt.show()