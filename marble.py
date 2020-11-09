# marble.py

import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from scipy.misc import derivative
import argparse


################################### START ARGUMENTS ###################################

parser = argparse.ArgumentParser(description='Script purpose: Simulate and plot a marble in a hole')

parser.add_argument('--mass', type=int, default=1, required=False,
                     help='mass. in units of g, default is 1g.')

parser.add_argument('--start', type=list, default=np.array([0,0]), required=False,
                     help='start. the starting point of the marble. default is ([0,0]).')

parser.add_argument('--v0', type=np.array, default=np.array([0,0]), required=False,
                     help='v0. Initial speed of the marble in format of ([0,0]). units of m/s, default is 0m/s.')

parser.add_argument('--intervals', type=int, default=0.1, required=False,
                     help='intervals. the time gaps between snapshots of the plot. units of s, default is 100ms.')

parser.add_argument('--points', type=int, default=50, required=False,
                     help='points. the number of snapshots taken, default is 50 points.')

args = parser.parse_args()

mass = args.mass
start = args.start
v0 = args.v0
intervals = args.intervals
points = args.points

################################### END ARGUMENTS ###################################

################################### CLASS FUNCTIONS #################################

class marble:
    def __init__(self, mass, start, v0, intervals, points, *args, **kwargs):
        self.mass = mass
        self.loc = start
        self.speed = v0
        self.intervals = intervals
        self.intervals_count = 0
        self.points = points
        self.inarena = 1


    #creating the shape of the hole.
    def fsin(self,x):
        return -np.sin(x)

    def create_arena(self):
        x=np.linspace(0,np.pi)
        return plt.plot(x,marble.fsin(x))

    def get_slope(self,x):
        deriv = derivative(marble.fsin, x, dx=1e-9)
        return deriv

    # calculate the force vector of the gravity.
    def get_Fg(self):
        return np.array([self.loc[0],-self.mass*9.8])


    # Newton's equations
    def get_A(self, force):
        # adjust direction of A
        force[0] = np.linalg.norm(force) * -marble.get_slope(self.loc[0])
        force[1] = marble.fsin(force[0])
        return force/(self.mass)

    def get_V(self, A):
        return self.speed + A * self.intervals

    def get_X(self, V, A):
        new_X = self.loc + V * self.intervals + 0.5 * A * pow((self.intervals),2)
        if not marble.is_next_valid(new_X):
            self.inarena = 0
        new_X[1] = marble.fsin(new_X[0])      
        return new_X


    # check if the marble is still in the arena
    def is_next_valid(self, X):
        if X[0] > np.pi or X[0] < 0 :
            return False
        return True

    # main simulation function
    def simulate(self):

        marble.create_arena()

        while self.intervals_count < self.points:
            new_A = marble.get_A(marble.get_Fg())

            new_V = marble.get_V(new_A)

            new_X = marble.get_X(new_V, new_A)

            # case the marble flew out
            if self.inarena == 0:
                plt.text(0, -1.2, "mass flew out of arena after "+ str(self.intervals * self.intervals_count) + " sec.")
                plt.show()
                return

            add_data = plt.plot(new_X[0], new_X[1], 'r+')
            plt.text(new_X[0],new_X[1], str(self.intervals_count))
            self.speed = new_V
            self.loc = new_X
            self.intervals_count += 1

        plt.show()
        return
            

marble = marble(mass, start, v0, intervals, points)
marble.simulate()


        