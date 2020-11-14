# marble.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
import argparse
import progressbar
from datetime import datetime
import os

################################### START ARGUMENTS ###################################

parser = argparse.ArgumentParser(description='Script purpose: Simulate and plot the position of a marble in a hole')

parser.add_argument('--mass', type=int, default=1, required=False,
                     help='mass of the marble. units of g, default is 1g.')

parser.add_argument('--start', type=str, default='0,0', required=False,
                     help='Starting point of the marble. format: start=x,y , default is [0,0].')

parser.add_argument('--v0', type=str, default='0,0', required=False,
                     help='Initial speed vector of the marble. format: v0=x,y , default is [0,0]m/s.')

parser.add_argument('--intervals', type=int, default=0.1, required=False,
                     help='Time gaps between snapshots of the plot. units of s, default is 0.1s.')

parser.add_argument('--points', type=int, default=100, required=False,
                     help='Number of snapshots taken, default is 100 points.')

parser.add_argument('--save_copy', type=str, default='True', required=False,
                     help='Saves a copy of the output. default is False.')

parser.add_argument('--show', type=str, default='True', required=False,
                     help='Show me the results. default is True. If false, saves a copy of the results.')

args = parser.parse_args()

mass = args.mass
start = np.fromstring(args.start, sep=',')
v0 = np.fromstring(args.v0, sep=',')
intervals = args.intervals
points = args.points
save_copy = args.save_copy
if save_copy.lower() in ['false','f']:
    save_copy = False
else:
    save_copy = True
show = args.show
if show.lower() in ['false','f']:
    show = False
else:
    show = True

if not show:
    save_copy = True
    print("\n\tNOTE: simulation output will be saved since you chose not to show results.")

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
        plt.text(-0.1, -0.7, "points: " + str(self.intervals_count-1))
        plt.text(-0.1, -0.8, "intervals: " + str(self.intervals)+"s")
        plt.text(-0.1, -0.9, "duration: " + str(format(self.intervals * (self.intervals_count-1),'.3f')) + "s")
        plt.text(-0.1, -1, "v0: (" +str(format(v0[0],'.3f'))+","+str(format(v0[1],'.3f'))+") m/s")


    def save_plots(self, name):
        if not os.path.exists('output_plots'):
            os.makedirs('output_plots')
        path = 'output_plots\\'+str(datetime.now().strftime("%x")).replace('/','-')+'_'+ str(datetime.now().strftime("%X")).replace(':','')+'_'+str(name)+'.png'
        plt.savefig(path)
        print("\n\tSaved "+name+" plot at:\n\t"+path)

    # main simulation function
    def simulate(self):

        print("\n\tSetting up marble simulation")
        print("\tmass: " + str(self.mass) + "g")
        print("\tintervals: " + str(self.intervals)+"s")
        print("\tv0: (" +str(format(v0[0],'.3f'))+" , "+str(format(v0[1],'.3f'))+") m/s")
        print("\tstating point: (" +str(format(start[0],'.3f'))+" , "+str(format(start[1],'.3f'))+")\n")

        # lists to collect data for X(t)
        t_list =[]
        x_list =[]
        t_list.append(0)
        x_list.append(start[0])

        marble.create_arena()

        # setteing up Y(t) plot
        plot_XY = plt.figure(1)
        plt.suptitle('Marble simulation output: Y(x)', fontsize=14, fontweight='bold')
        plt.xlabel("X [m]")
        plt.ylabel("Y [m]")
        plt.text(self.loc[0],self.loc[1], "_______START("+str(format(self.loc[0],'.3f'))+","+str(format(self.loc[1],'.3f'))+")", color='green', fontsize=12)
        add_data = plt.plot(self.loc[0], self.loc[1], 'gX')

        print("\tarena created.\n")
        starttime = datetime.now()
        print("\tSimulation in progress...")

        # creating progress bar
        bar = progressbar.ProgressBar(maxval=self.points, \
        widgets=[progressbar.Bar('=', '\t[', ']'), ' ', progressbar.Percentage()])
        bar.start()

        # main simulation loop
        while self.intervals_count <= self.points:

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
            bar.update(self.intervals_count)
            self.speed = new_V
            self.loc = new_X
            self.intervals_count += 1
            t_list.append(self.intervals_count*self.intervals)
            x_list.append(self.loc[0])
        
        bar.finish()
        marble.add_data_to_plot()
        print("\n\tSimulation finished after "+str(datetime.now()-starttime)+"\n")

        if self.intervals_count > 2000:
            print("\n\tPreparing data, please wait..."+"\n")
        if save_copy:
            marble.save_plots('Yx')

        # build plot X(t)
        plot_XT= plt.figure(2)
        plt.scatter(t_list, x_list)
        plt.plot(t_list, x_list)
        plt.suptitle('Marble simulation output: X(t)', fontsize=14, fontweight='bold')
        plt.xlabel("t [s]")
        plt.ylabel("X [m]")
        plt.text(t_list[0], x_list[0], "    START\n    (x="+str(format(x_list[0],'.3f'))+"m)", color='green', fontsize=10)
        plt.text(t_list[-1], x_list[-1], "END\n(x="+str(format(x_list[-1],'.3f'))+"m)", color='red', fontsize=10)

        if save_copy:
            marble.save_plots('Xt')
        else:
            print ("\n\tNOTE: you can autosave the plots next time by adding 'save_copy=t' to the simulation cmd.\n")

        if show:
            print ("\n\tDisplaying plot Y(x):\n")
            print ("\tDisplaying plot X(t):\n")
            print ("\tClose all plot windows and press ENTER to continue...\n")
        print("\n\tScript finished after "+str(datetime.now()-starttime)+"\n")
        plt.show(block=show)
        print ("Done.")
        return
            
################################### END FUNCTIONS ###################################

#################################### EXECUTION ######################################

marble = marble(mass, start, v0, intervals, points)
marble.simulate()


        