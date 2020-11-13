# marble.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
import argparse
from datetime import datetime


################################### START ARGUMENTS ###################################

parser = argparse.ArgumentParser(description='Script purpose: Simulate and plot a marble in a hole')

parser.add_argument('--mass', type=int, default=1, required=False,
                     help='mass. in units of g, default is 1g.')

parser.add_argument('--start', type=list, default=np.array([0,0]), required=False,
                     help='start. the starting point of the marble. default is ([0,0]).')

parser.add_argument('--v0', type=np.array, default=np.array([0,0]), required=False,
                     help='v0. Initial speed of the marble in format of ([0,0]). Vector size in units of m/s, default is 0m/s.')

parser.add_argument('--intervals', type=int, default=0.1, required=False,
                     help='intervals. the time gaps between snapshots of the plot. units of s, default is 100ms.')

parser.add_argument('--points', type=int, default=100, required=False,
                     help='points. the number of snapshots taken, default is 100 points.')

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
        return np.array([self.loc[0], self.loc[1]-self.mass*9.8])


    def get_Fu(self):
        liquid_viscosity = 0.89
        shear_rate = 0.00001*(1-abs(self.loc[1]))
        area_of_plate = abs(np.pi / 2 -self.loc[0]) * 2
        Fu = liquid_viscosity * shear_rate * area_of_plate
        ret=np.array([self.loc[0] - abs(self.speed[0] * Fu) / self.speed[0] , self.loc[1] - abs(self.speed[1] * Fu) / self.speed[1]])
        return ret

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

    # add frame and parameter data to plot
    def add_data_to_plot(self):
        plt.text(-0.1, -0.6, "mass: " + str(self.mass) + "g")
        plt.text(-0.1, -0.7, "points: " + str(self.intervals_count))
        plt.text(-0.1, -0.8, "intervals: " + str(self.intervals)+"s")
        plt.text(-0.1, -0.9, "duration: " + str(format(self.intervals * (self.intervals_count),'.3f')) + "s")
        plt.text(-0.1, -1, "v0: (" +str(format(v0[0],'.3f'))+","+str(format(v0[1],'.3f'))+") m/s")

    # main simulation function
    def simulate(self):

        print("Setting up marble simulation")
        print("mass: " + str(self.mass) + "g")
        print("intervals: " + str(self.intervals)+"s")
        print("v0: (" +str(format(v0[0],'.3f'))+","+str(format(v0[1],'.3f'))+") m/s")
        print("stating point: (" +str(format(start[0],'.3f'))+","+str(format(start[1],'.3f'))+")\n")

        # lists to collect data for X(t)
        t_list =[]
        x_list =[]
        t_list.append(0)
        x_list.append(start[0])

        marble.create_arena()
        plot_XY = plt.figure(1)
        plt.suptitle('Marble simulation output: Y(x)', fontsize=14, fontweight='bold')
        plt.xlabel("X [m]")
        plt.ylabel("Y [m]")
        plt.text(self.loc[0],self.loc[1], "_______START("+str(format(self.loc[0],'.3f'))+","+str(format(self.loc[1],'.3f'))+")", color='green', fontsize=12)
        add_data = plt.plot(self.loc[0], self.loc[1], 'gX')

        print("arena created.\n")
        starttime = datetime.now()
        print("Starting simulation...")
        while self.intervals_count <= self.points:

            Feq = np.add(marble.get_Fg(), marble.get_Fu())

            new_A = marble.get_A(marble.get_Fg())

            new_V = marble.get_V(new_A)

            new_X = marble.get_X(new_V, new_A)

            # case the marble flew out
            if self.inarena == 0:
                plt.text(0, 0.1, "mass flew out of the arena after "+ str(self.intervals * self.intervals_count) + " sec.")
                break
            
            if self.intervals_count == self.points:
                add_data = plt.plot(new_X[0], new_X[1], 'rX')
                plt.text(new_X[0],new_X[1], "________END("+str(format(new_X[0],'.3f'))+","+str(format(new_X[1],'.3f'))+")", color='red', fontsize=12)

                #plt.text(0, 0.1, "Simulation ran for "+ str(self.intervals * self.intervals_count) + " sec, ploting "+ str(self.points) + " points at "+str(self.intervals)+"s intervals.")
            else:
                add_data = plt.plot(new_X[0], new_X[1], 'g+')
                plt.text(new_X[0],new_X[1], str(self.intervals_count + 1))
            self.speed = new_V
            self.loc = new_X
            self.intervals_count += 1
            t_list.append(self.intervals_count*self.intervals)
            x_list.append(self.loc[0])
        
        marble.add_data_to_plot()
        print("Finished after "+str(datetime.now()-starttime)+"\n")

        # build plot X(t)
        plot_XT= plt.figure(2)
        plt.scatter(t_list, x_list)
        plt.plot(t_list, x_list)
        plt.suptitle('Marble simulation output: X(t)', fontsize=14, fontweight='bold')
        plt.xlabel("t [s]")
        plt.ylabel("X [m]")
        plt.text(t_list[0], x_list[0], "    START\n    (x="+str(format(x_list[0],'.3f'))+"m)", color='green', fontsize=10)
        plt.text(t_list[-1], x_list[-1], "END\n(x="+str(format(x_list[-1],'.3f'))+"m)", color='red', fontsize=10)

        print ("Displaying plot Y(x):\n")
        print ("Displaying plot X(t):\n")
        plt.show()
        print ("Done.")
        return
            
################################### END FUNCTIONS ###################################

#################################### EXECUTION ######################################

marble = marble(mass, start, v0, intervals, points)
marble.simulate()


        