import emit.emitting_particles as ep
import emit.search_sigmas as search_sigmas
import numpy as np
import math
import matplotlib.pyplot as plt

# mass of electron
m = 9.1 * pow(10, -31)
# speed of light
c = 299792458
# eV in J
eV = 1.6 * pow(10, -19)
# Mega
Mega = 1000000
# min energy in MeV
E_min = 0.001
# min virtual mass of particle
W_min = pow(10, -12)
# sigma Compton
sigma_c = 0.3107
# sigma of rest interferences
sigma_r = 0.1624
# probability of Compton scattering
prob_c = sigma_c / ep.sigma
# param of chemical element
Z = 82
r_e = 2.81794 * pow(10, -13)
# number of detectors at one axis
n = 100
# central detectors
d_central = np.linspace(0, ep.h, n)
density_central = [0] * n
# Sides detectors
z_sides = np.linspace(0, ep.h, n)
x_left_side = ep.r
x_right_side = -ep.r
density_left_side = [0] * n
density_right_side = [0] * n
# sides at the top
z_top_sides = ep.h
x_top_a_side = np.linspace(-ep.r, ep.r, n)
y_top_b_side = np.linspace(-ep.r, ep.r, n)
density_top_a_side = [0] * n
density_top_b_side = [0] * n

# initial parameters ( energy in MeV)
E_0 = 3
W_0 = 1

ep.emitter(False)


# ep.drawing()


# Get initial guide cos
def getGuideCos(x1, y1, z1, x0, y0, z0):
    l = np.sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2) + pow(z1 - z0, 2))
    return np.abs(x1 - x0) / l, np.abs(y1 - y0) / l, np.abs(z1 - z0) / l


# Get point after interference
def getFinalPoint(x1, y1, z1, x0, y0, z0, cos_teta, sigma):
    init_cos = getGuideCos(x1, y1, z1, x0, y0, z0)
    phi = np.random.uniform(0, 2 * np.pi)
    w3 = init_cos[2] * cos_teta + np.sqrt((1 - pow(cos_teta, 2)) * (1 - pow(init_cos[2], 2))) * np.cos(phi)
    w2 = (init_cos[1] * (cos_teta - w3 * init_cos[2]) + init_cos[0] * np.sin(phi) * np.sqrt(
        (1 - pow(cos_teta, 2)) * (1 - pow(init_cos[2], 2)))) / (1 - pow(init_cos[2], 2))
    w1 = (init_cos[0] * (cos_teta - w3 * init_cos[2]) - init_cos[1] * np.sin(phi) * np.sqrt(
        (1 - pow(cos_teta, 2)) * (1 - pow(init_cos[2], 2)))) / (1 - pow(init_cos[2], 2))
    l = -1 / sigma * math.log(np.random.random())
    return l * w1, l * w2, l * w3


for i in range(ep.n):
    x0 = ep.x0[i]
    y0 = ep.y0[i]
    z0 = 0
    x1 = ep.x1[i]
    y1 = ep.y1[i]
    z1 = ep.z1[i]
    E = E_0
    W = W_0
    while ep.isInArea(x1, y1, z1) and np.random.random() <= prob_c:
        teta = np.random.uniform(0, np.pi / 2)
        cos_teta = np.cos(teta)
        alpha = E / (m * pow(c, 2))
        E = E * pow((1 + E * Mega * eV * (1 - cos_teta) / (m * pow(c, 2))), -1)
        sigmas = search_sigmas.findSigma(E)
        W = W * sigmas[0] / sigmas[1]
        if W < W_min and E < E_min:
            continue
        sigma_int = 2 * np.pi * Z * pow(r_e, 2) * (
            (1 + alpha) / pow(alpha, 2) * (2 * (1 + alpha) / (1 + 2 * alpha) - np.log(1 + 2 * alpha) / alpha) + np.log(
                1 + 2 * alpha) / (2 * alpha) - (1 + 3 * alpha) / pow(1 + 2 * alpha, 2))
        sigma_dif = Z * pow(r_e, 2) / 2 * pow(1 + alpha * (1 - cos_teta), -2) * (
            1 + pow(cos_teta, 2) + (pow(alpha, 2) * pow(1 - cos_teta, 2)) / (1 + alpha * (1 - cos_teta)))
        k = sigma_dif / sigma_int
        for j in range(n):
            R_central = np.sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1 - d_central[j], 2))
            density_central[j] = density_central[j] + k * W * np.exp(-sigmas[1] * R_central) / pow(R_central, 2)
            R_left_side = np.sqrt(pow(x1 - x_left_side, 2) + pow(y1, 2) + pow(z1 - z_sides[j], 2))
            R_right_side = np.sqrt(pow(x1 - x_right_side, 2) + pow(y1, 2) + pow(z1 - z_sides[j], 2))
            density_left_side[j] = density_left_side[j] + k * W * np.exp(-sigmas[1] * R_left_side) / pow(R_left_side, 2)
            density_right_side[j] = density_right_side[j] + k * W * np.exp(-sigmas[1] * R_right_side) / pow(
                R_right_side, 2)
            R_top_a_side = np.sqrt(pow(x1 - x_top_a_side[j], 2) + pow(y1, 2) + pow(z1 - z_top_sides, 2))
            R_top_b_side = np.sqrt(pow(x1, 2) + pow(y1 - y_top_b_side[j], 2) + pow(z1 - z_top_sides, 2))
            density_top_a_side[j] = density_top_a_side[j] + k * W * np.exp(-sigmas[1] * R_top_a_side) / pow(
                R_top_a_side, 2)
            density_top_b_side[j] = density_top_b_side[j] + k * W * np.exp(-sigmas[1] * R_top_b_side) / pow(
                R_top_b_side, 2)
        final_point = getFinalPoint(x1, y1, z1, x0, y0, z0, cos_teta, sigmas[1])
        x0 = x1
        y0 = y1
        z0 = z1
        x1 = final_point[0] + x1
        y1 = final_point[1] + y1
        z1 = final_point[2] + z1

        # ep.ax.plot([x0, x1], [y0, y1], [z0, z1], linewidth=1, color='red')


def drawing():
    fig_central = plt.figure()
    central_plot = fig_central.add_subplot(111)
    central_plot.set_title("Ось симметрии")
    fig_left_side = plt.figure()
    left_side_plot = fig_left_side.add_subplot(111)
    left_side_plot.set_title("Левая стенка")
    fig_right_side = plt.figure()
    right_side_plot = fig_right_side.add_subplot(111)
    right_side_plot.set_title("Правая стенка")
    fig_top_a_side = plt.figure()
    top_a_side_plot = fig_top_a_side.add_subplot(111)
    top_a_side_plot.set_title("Основание а")
    fig_top_b_side = plt.figure()
    top_b_side_plot = fig_top_b_side.add_subplot(111)
    top_b_side_plot.set_title("Основание b")
    central_plot.plot(d_central, density_central)
    right_side_plot.plot(z_sides, density_right_side)
    left_side_plot.plot(z_sides, density_left_side)
    top_a_side_plot.plot(x_top_a_side, density_top_a_side)
    top_b_side_plot.plot(y_top_b_side, density_top_b_side)
    plt.show()


drawing()
