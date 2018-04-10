import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

# Parameters of equation
a = 1.5
l = 2 * np.pi
T = 2
# Steps
N = 1000
M = 100
x = np.linspace(0, l, N)
t = np.linspace(0, T, M)
h = x[1]
tau = t[1]
# Condition of pseudo stability
print(tau <= pow(h, 2) / (2 * pow(a, 2)))
# Temperature
u = []
# Initial condition
u.append([3 * (l - i) for i in x])


# Function
def f(x, t):
    return t * np.sin(x)


'''Implicit method'''
a_i = tau * pow(a, 2) / pow(h, 2)
b_i = -(1 + 2 * tau * pow(a, 2) / pow(h, 2))
c_i = a_i
for j in range(M - 1):
    A = [1]
    B = [0]
    for i in range(1, N):
        d = -tau*f(x[i], t[j + 1]) - u[j][i]
        A.append(-c_i / (b_i + a_i * A[i - 1]))
        B.append((d - a_i * B[i - 1]) / (b_i + a_i * A[i - 1]))
    y = [0] * N
    y[N - 1] = B[N - 1] / (1 - A[N - 1])
    for i in reversed(range(N - 1)):
        y[i] = A[i] * y[i + 1] + B[i]
    u.append(y)

fig, ax = plt.subplots()

line, = ax.plot(x, u[0])

def animate(i):
    line.set_ydata(u[i])  # update the data
    return line,


ani = animation.FuncAnimation(fig, animate, interval=100, frames=M, repeat=False)

plt.show()