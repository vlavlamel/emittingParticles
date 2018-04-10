import matplotlib.pyplot as plt
import numpy as np

eps = 0.01
r = 30

# Testing binary method
'''x_circle = np.linspace(-r, r, 100)
y_circle = [np.sqrt(pow(r, 2) - pow(i, 2)) for i in x_circle]
fig1 = plt.figure()
ax = fig1.add_subplot(111)'''


def getCrossoverPoint(x0, y0, x1, y1):
    global eps, x_circle, y_circle
    x = 0
    y = 0
    start_x = x0
    start_y = y0
    end_x = x1
    end_y = y1
    while abs(pow(x, 2) + pow(y, 2) - pow(r, 2)) > eps:
        x = (end_x - start_x) / 2 + start_x
        y = (end_y - start_y) / 2 + start_y
        if pow(x, 2) + pow(y, 2) > pow(r, 2):
            end_x = x
            end_y = y
        else:
            start_x = x
            start_y = y
    '''ax.plot(x0, y0, '-o', x1, y1, '-o', x, y, '-o')
    ax.plot((x0, x1), (y0, y1), '-b')
    ax.plot(x_circle, y_circle)
    plt.show()'''
    # print(np.sqrt(pow(x, 2) + pow(y, 2)))
    return x, y, np.sqrt(pow(x - x0, 2) + pow(y - y0, 2))


# getCrossoverPoint(20, 20, -35, 0)