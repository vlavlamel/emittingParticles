import matplotlib.pyplot as plt
import numpy as np
import math
from emit.search_crossover_point import getCrossoverPoint
from mpl_toolkits.mplot3d import Axes3D

# fig1 = plt.figure()
# fig2 = plt.figure()
# ax = fig1.add_subplot(111, projection='3d')
# ad = fig2.add_subplot(111)

# Parameters of cylinder
r = 30
h = 20
# Parameters of ellipse
a = 20
b = 10
# Parameter of mean free path

# without density
#sigma = 0.04172

# with density
sigma = 0.4731

# Quantity of particles
n = 1000
# List of initial points
x0 = []
y0 = []
# List of final points
x1 = []
y1 = []
z1 = []


# Get initial points
def initialPoints():
    global n, x0, y0, b, a
    counter = 0
    while counter < n:
        y = np.random.uniform(0, 2 * a)
        x = np.random.uniform(-a, a)
        if (y <= b * np.sqrt(1 - x ** 2 / a ** 2)):
            counter += 1
            x0.append(x)
            y0.append(y)


# Get final points
def emitter(setInBound):
    global h, n, r, x1, y1, z1
    for i in range(n):
        phi = 2 * np.pi * np.random.random()
        cosTeta = np.random.random()
        l = -1 / sigma * math.log(np.random.random())
        x = x0[i] + l * np.sqrt(1 - pow(cosTeta, 2)) * np.cos(phi)
        y = y0[i] + l * np.sqrt(1 - pow(cosTeta, 2)) * np.sin(phi)
        z = l * cosTeta
        if setInBound:
            if z > h:
                Lxy = np.sqrt(1 - pow(cosTeta, 2)) / cosTeta * h
                x = x0[i] + Lxy * np.cos(phi)
                y = y0[i] + Lxy * np.sin(phi)
                z = h
            if pow(x, 2) + pow(y, 2) > pow(r, 2):
                boundary = getCrossoverPoint(x0[i], y0[i], x, y)
                x = boundary[0]
                y = boundary[1]
                z = boundary[2] * cosTeta / np.sqrt(1 - pow(cosTeta, 2))
        x1.append(x)
        y1.append(y)
        z1.append(z)


# Draw hedgehog
def drawing():
    global a, b, r, h, n, x0, y0, x1, y1, z1  # , ax, ad
    fig1 = plt.figure()
    fig2 = plt.figure()
    ax = fig1.add_subplot(111, projection='3d')
    ad = fig2.add_subplot(111)
    #Ellipse
    x_ellipse = np.linspace(-a, a, 100)
    y_ellipse = b * np.sqrt(1 - x_ellipse ** 2 / a ** 2)
    ad.plot(x_ellipse, y_ellipse, linewidth=1, color='blue')
    ad.plot(x_ellipse, np.zeros(x_ellipse.size), linewidth=1, color='blue')
    ad.plot(x0, y0, '.')

    # Cylinder
    x_cylinder = np.linspace(-r, r, 100)
    z_cylinder = np.linspace(0, h, 100)
    Xc, Zc = np.meshgrid(x_cylinder, z_cylinder)
    Yc = r * np.sqrt(1 - (Xc / r) ** 2)

    # Lines
    for i in range(n):
        X = np.array([x0[i], x1[i]])
        Y = np.array([y0[i], y1[i]])
        Z = np.array([0, z1[i]])
        ax.plot(X, Y, Z, linewidth=1, color='green')

    # Draw parameters
    rstride = 20
    cstride = 10
    ax.plot_surface(Xc, Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride, color='blue')
    ax.plot_surface(Xc, -Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride, color='blue')
    ax.plot(x_ellipse, y_ellipse, 0, linewidth=1, color='blue')
    ax.plot(x_ellipse, np.zeros(x_ellipse.size), 0, linewidth=1, color='blue')

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.show()


# Return boolean - is particle in cylinder
def isInArea(x, y, z):
    global r, h
    if (z > h) or (z < 0) or (pow(x, 2) + pow(y, 2) > pow(r, 2)):
        return False
    else:
        return True


initialPoints()
#emitter(True)
#drawing()