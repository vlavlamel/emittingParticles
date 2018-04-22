import numpy as np
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# Parameters of equation
Do = 0.01
k = 0.05
beta = 5 * pow(10, -7)
thickness = 100
H1 = pow(10, 5)
H2 = 2 * pow(10, 5)
T = 1000
# Steps
N = 1000
M = 1000
x = np.linspace(0, thickness, M)
t = np.linspace(0, T, N)
h = x[1]
tau = t[1]
# Concentration
u = []


# Diffusion coefficient
def D(x):
    global Do, thickness
    return Do * (1 + 2 * x / thickness)


#
def a(i):
    global x
    return D((x[i] - x[i - 1]) / 2)


def a_next(i):
    global x
    return D((x[i + 1] - x[i]) / 2)


# Tridiagonal matrix algorithm
def tridiagonal():
    u.clear()
    # Initial condition
    u.append([1] * M)
    for j in range(N - 1):
        # A = [1]
        A = [-(a_next(0) / h) / (-k - a_next(0) / h)]
        B = [0]
        for i in range(1, M - 1):
            a_i = -tau * a(i) / pow(h, 2)
            b_i = 1 + tau * a_next(i) / pow(h, 2) + tau * a(i) / pow(h, 2) + tau * beta
            c_i = -tau * a_next(i) / pow(h, 2)
            d = u[j][i]
            A.append(-c_i / (b_i + a_i * A[i - 1]))
            B.append((d - a_i * B[i - 1]) / (b_i + a_i * A[i - 1]))
        y = [0] * M
        Am = a(M - 1) / (a(M - 1) + k * h)
        y[M - 1] = Am * B[M - 2] / (1 - Am * A[M - 2])
        for i in reversed(range(M - 1)):
            y[i] = A[i] * y[i + 1] + B[i]
        u.append(y)

tridiagonal()
'''Animation'''
fig, ax = plt.subplots()

line, = ax.plot(x, u[0])
ax.set_yticks(np.arange(0, 1.1, 0.1))


def animate(i):
    line.set_ydata(u[i])  # update the data
    return line,


ani = animation.FuncAnimation(fig, animate, interval=1, frames=N, repeat=False)
plt.show()

'''3d plot'''
fig1 = plt.figure()
ac = fig1.gca(projection='3d')

X, Y = np.meshgrid(x, t)
surf = ac.plot_surface(X, Y, u, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ac.set_xlabel('X')
ac.set_ylabel('t')
ac.set_zlabel('C')

# tridiagonal()
# plt.figure(1)
# plt.plot(x, u[N - 1])
# plt.yticks(np.arange(0, 1.1, 0.1))
# plt.title("Do = 0.01, beta = 5 10^(-7), k = 0.05")
#
# Do = Do * 5
# tridiagonal()
# plt.figure(2)
# plt.plot(x, u[N - 1])
# plt.yticks(np.arange(0, 1.1, 0.1))
# plt.title("Do = 0.05, beta = 5 10^(-7), k = 0.05")
#
# Do = Do / 5
# beta = beta * 1000
# tridiagonal()
# plt.figure(3)
# plt.plot(x, u[N - 1])
# plt.yticks(np.arange(0, 1.1, 0.1))
# plt.title("Do = 0.01, beta = 5 10^(-4), k = 0.05")
#
# beta = beta / 1000
# k = k / 100
# tridiagonal()
# plt.figure(4)
# plt.plot(x, u[N - 1])
# plt.yticks(np.arange(0, 1.11, 0.1))
# plt.title("Do = 0.01, beta = 5 10^(-7), k = 0.0005")

plt.show()
