import numpy as np
import matplotlib.pyplot as plt
import math
import time

SIMULATION_TIME = 0.005
k = 9e9
G = 6.67e-11


def calculate_movement(x, v, a, t):
    return x + v * t + 0.5 * a * t ** 2


def calculate_velocity(v, a, t):
    return v + a * t


def cartesian_distance(x1, y1, x2, y2):
    return ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5


class Charge:
    def __init__(self, m, x, y, c, vx, vy, ax, ay):
        """
        :param m: mass
        :param x: x position
        :param y: y position
        :param c: charge
        :param vx: velocity on x axis
        :param vy: velocity on y axis
        :param ax: acceleration on x axis
        :param ay: acceleration on y axis
        """
        self._m = m
        self._x = x
        self._y = y
        self._c = c
        self._vx = vx
        self._vy = vy
        self._ax = ax
        self._ay = ay

    def update_location(self, t):
        self._x = calculate_movement(self._x, self._vx, self._ax, t)
        self._y = calculate_movement(self._y, self._vy, self._ay, t)

    def update_velocity(self, t):
        self._vx = calculate_velocity(self._vx, self._ax, t)
        self._vy = calculate_velocity(self._vy, self._ay, t)

    def update_acceleration(self, ax, ay):
        """
        Updates the acceleration to the given values.
        :param ax: The acceleration on the x axis.
        :param ay: The acceleration on the y axis.
        """
        self._ax = ax
        self._ay = ay

    def update_acceleration_by_electric_force(self, other):
        # Calculate the size of the force
        r = cartesian_distance(self._x, self._y, other._x, other._y)
        if r == 0:
            return

        f = k * self._c * other._c / r ** 2

        # Calculate the direction of the force
        if self._x == other._x or self._y == other._y:
            a = math.pi
        else:
            a = math.atan(abs(self._y - other._y) / abs(self._x - other._x))

        if self._x >= other._x and self._y <= other._y:
            a = -a
        elif self._x <= other._x and self._y <= other._y:
            a = math.pi + a
        elif self._x <= other._x and self._y >= other._y:
            a = math.pi - a

        # Calculate the force of each axis
        fx = f * math.cos(a)
        fy = f * math.sin(a)

        # Update the accelerations
        self._ax += fx / self._m
        self._ay += fy / self._m

    def update_acceleration_by_gravitation_force(self, other):
        # Calculate the size of the force
        r = cartesian_distance(self._x, self._y, other._x, other._y)
        if r == 0:
            return

        f = G * self._m * other._m / r ** 2

        # Calculate the direction of the force
        if self._x == other._x or self._y == other._y:
            a = math.pi
        else:
            a = math.atan(abs(self._y - other._y) / abs(self._x - other._x))

        if self._x <= other._x and self._y >= other._y:
            a = -a
        elif self._x >= other._x and self._y >= other._y:
            a = math.pi + a
        elif self._x >= other._x and self._y <= other._y:
            a = math.pi - a

        # Calculate the force of each axis
        fx = f * math.cos(a)
        fy = f * math.sin(a)

        # Update the accelerations
        self._ax += fx / self._m
        self._ay += fy / self._m

    def update(self, t):
        self.update_location(t)
        self.update_velocity(t)

def plot(charges, scatters):
    for i, charge in enumerate(charges):
        # if len(scatters) > i:
        scatters[i].remove()
        scatters[i] = plt.scatter(charge._x, charge._y, color='red' if charge._c < 0 else 'blue')

    plt.pause(SIMULATION_TIME)


def update(charges):
    for charge in charges:
        charge.update_acceleration(0, 0)
        for other in charges:
            if other is charge:
                continue
            charge.update_acceleration_by_electric_force(other)
            # charge.update_acceleration_by_gravitation_force(other)
        charge.update(SIMULATION_TIME)


def main():
    # Generate charges
    charges = []
    charges.append(Charge(0.05, 0.125, 0, -2e-6, 0, 0, 0, 0))
    charges.append(Charge(0.05, -0.125, 0, 8e-6, 0, 0, 0, 0))
    charges.append(Charge(0.05, 0, 0.2165, 3e-6, 0, 0, 0, 0))
    charges.append(Charge(0.05, 1, 0.4, 4e-6, 0, 0, 0, 0))
    charges.append(Charge(0.05, 1.5, -0.2, -7e-6, 0, 0, 0, 0))

    # Generate scatters list
    scatters = []
    for charge in charges:
        scatters.append(plt.scatter(charge._x, charge._y))

    plt.axis([-3, 3, -3, 3])
    # Maximise the plotting window
    # plot_backend = plt.get_backend()
    # mng = plt.get_current_fig_manager()
    # if plot_backend == 'TkAgg':
    #     mng.window.wm_geometry("+0+0")
    #     mng.resize(*mng.window.maxsize())
    # elif plot_backend == 'wxAgg':
    #     mng.window.SetPosition((500, 0))
    #     mng.frame.Maximize(True)
    # elif plot_backend == 'Qt4Agg':
    #     mng.window.showMaximized()


    while (True):
        update(charges)
        plot(charges, scatters)


if __name__ == '__main__':
    main()
