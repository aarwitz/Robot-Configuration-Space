"""
Classes to define potential and potential planner for the sphere world
"""

import numpy as np
import me570_geometry
from matplotlib import pyplot as plt
from scipy import io as scio
import math

class SphereWorld:
    """ Class for loading and plotting a 2-D sphereworld. """
    def __init__(self):
        """
        Load the sphere world from the provided file sphereworld.mat, and sets the
    following attributes:
     -  world: a  nb_spheres list of  Sphere objects defining all the spherical obstacles in the
    sphere world.
     -  x_start, a [2 x nb_start] array of initial starting locations (one for each column).
     -  x_goal, a [2 x nb_goal] vector containing the coordinates of different goal locations (one
    for each column).
        """
        data = scio.loadmat('sphereWorld.mat')

        self.world = []
        for sphere_args in np.reshape(data['world'], (-1, )):
            sphere_args[1] = np.asscalar(sphere_args[1])
            sphere_args[2] = np.asscalar(sphere_args[2])
            self.world.append(me570_geometry.Sphere(*sphere_args))

        self.x_goal = data['xGoal']
        self.x_start = data['xStart']
        self.theta_start = data['thetaStart']

    def plot(self):
        """
        Uses Sphere.plot to draw the spherical obstacles together with a  * marker at the goal location.
        """

        for sphere in self.world:
            sphere.plot('r')
        ax = plt.gca()
        ax.scatter(self.x_goal[0, :], self.x_goal[1, :], c='g', marker='*')

        plt.xlim([-11, 11])
        plt.ylim([-11, 11])

class RepulsiveSphere:
    """ Repulsive potential for a sphere """
    def __init__(self, sphere):
        """
        Save the arguments to internal attributes
        """
        self.sphere = sphere

    def eval(self, x_eval):
        """
        Evaluate the repulsive potential from  sphere at the location x= x_eval. The function returns
    the repulsive potential as given by      (  eq:repulsive  ).
        """
        # sign_rad = np.sign(self.sphere.radius)
        
        sign_rad =1
        d_influence = self.sphere.distance_influence
        d_ix = self.sphere.distance(x_eval)[0]
        if d_ix > 0 and d_ix < d_influence:
            d_inverse = 1/d_ix - 1/d_influence
            u_rep = (1/2 * d_inverse**2)*sign_rad
        elif d_ix > d_influence:
            u_rep = 0
        else:
            u_rep = math.nan
        return u_rep

    def grad(self, x_eval):
        """
        Compute the gradient of U_rep for a single sphere, as given by (eq:repulsive-gradient).
        """
        d_influence = self.sphere.distance_influence
        d_ix = self.sphere.distance(x_eval)[0]
        # sign_rad = np.sign(self.sphere.radius)
        sign_rad=1
        if self.sphere.center[0] == 3.75:
            if d_ix > 0 and d_ix < d_influence:
                d_inverse = 1/d_ix - 1/d_influence
                grad_d_ix =self.sphere.distance_grad(self.sphere,x_eval)
                grad_u_rep = -1*(d_inverse) * 1/d_ix**2 * grad_d_ix*sign_rad
                # print('grad_u_rep')
                # print(grad_u_rep)
            elif d_ix > d_influence:
                grad_u_rep = np.zeros((2,1))
            else:
                grad_u_rep = math.nan
            return grad_u_rep
        else:
            if d_ix > 0 and d_ix < d_influence:
                d_inverse = 1/d_ix - 1/d_influence
                grad_d_ix =self.sphere.distance_grad(self.sphere,x_eval)
                grad_u_rep = -1*(d_inverse) * 1/d_ix**2 * grad_d_ix*sign_rad
            elif d_ix > d_influence:
                grad_u_rep = np.zeros((2,1))
            else:
                grad_u_rep = math.nan 
            return grad_u_rep

class Attractive:
    """ Repulsive potential for a sphere """
    def __init__(self, potential):
        """
        Save the arguments to internal attributes
        """
        self.potential = potential

    def eval(self, x_eval):
        """
        Evaluate the attractive potential  U_ attr at a point  xEval with respect to a goal location
    potential.xGoal given by the formula: If  potential.shape is equal to  'conic', use p=1. If
    potential.shape is equal to  'quadratic', use p=2.
        """
        x_eval = x_eval.reshape(2,1)
        x_goal = self.potential['x_goal']
        if self.potential['shape'] == 'conic':
            p = 1
        else:
            p = 2
        u_attr = np.linalg.norm(x_eval-x_goal)
        return u_attr**p

    def grad(self, x_eval):
        """
        Evaluate the gradient of the attractive potential  U_ attr at a point  xEval. The gradient is
    given by the formula If  potential['shape'] is equal to  'conic', use p=1; if it is equal to
    'quadratic', use p=2.
        """
        x_goal = self.potential['x_goal']
        vec_dist = x_eval-x_goal.reshape(2,1)
        dist = np.linalg.norm(vec_dist)
        u = self.eval(x_eval)
        if self.potential['shape'] == 'conic':
            p = 1
        else:
            p = 2
        # u = self.eval(x_eval)
        # grad_u_attr = p*u/u**(2)*vec_dist
        grad_u_attr = (p*u/dist**2)*vec_dist
        
        return grad_u_attr

class Total:
    """ Combines attractive and repulsive potentials """
    def __init__(self, world, potential):
        """
        Save the arguments to internal attributes
        """
        self.world = world
        self.potential = potential

    def eval(self, x_eval):
        """
        Compute the function U=U_attr+a*iU_rep,i, where a is given by the variable
    potential.repulsiveWeight
        """
        Attractive_eval = Attractive(self.potential).eval(x_eval)
        Repulsive_eval = sum([RepulsiveSphere(sphere).eval(x_eval) for sphere in self.world])
        u_eval = Attractive_eval + Repulsive_eval
        return u_eval

    def grad(self, x_eval):
        """
        Compute the gradient of the total potential,  U= U_ attr+    _i U_ rep,i, where   is given by
    the variable  potential.repulsiveWeight
        """
        Attractive_eval = Attractive(self.potential).grad(x_eval)
        alpha = self.potential['repulsiveWeight']
        Repulsive_eval = sum([alpha*RepulsiveSphere(sphere).grad(x_eval) for sphere in self.world])
        grad_u_eval = Attractive_eval + Repulsive_eval
        return grad_u_eval

class Planner:
    """  """
    def run(self, x_start, world, potential, planned_parameters):
        """
        This function uses a given function ( planner_parameters['control']) to implement a generic
    potential-based planner with step sierose  planner_parameters['epsilon'], and evaluates the cost
    along the returned path. The planner must stop when either the number of steps given by
    planner_parameters['nb_steps'] is reached, or when the norm of the vector given by
    planner_parameters['control'] is less than 5 10^-3 (equivalently,  5e-3).
        """
        x_start = x_start.reshape(2,1)
        max_nb_steps = planned_parameters['nb_steps']
        control = planned_parameters['control']
        epsilon = planned_parameters['epsilon']
        U = planned_parameters['U']
        x_path = x_start.reshape(2,1)
        u_path = np.array([U(x_start)])
        nb_steps = 0
        controlCurrent  = 1
        while nb_steps < max_nb_steps and np.linalg.norm(controlCurrent)>(5*10**(-3)):
        # while nb_steps < 10000 and u_path[-1]>(5*10**(-3)):
            nb_steps = nb_steps + 1
            controlCurrent = control(x_path[:,-1].reshape(2,1))
            x_updated = x_path[:,-1].reshape(2,1) - epsilon*controlCurrent.reshape(2,1)
            x_path = np.append(x_path,x_updated,axis=1)
            u = U(x_path[:,-1])
            u_path = np.append(u_path,u)
        # plt.show()
        while x_path.shape[1]<max_nb_steps:
            nan_col = np.array([[math.nan],[math.nan]])
            x_path = np.append(x_path,nan_col,axis=1)
            u_path = np.append(u_path,nan_col,axis=1)
        return x_path, u_path

    def run_plot(self):
        """
        This function performs the following steps:
     - Loads the problem data from the file !70!DarkSeaGreen2 sphereworld.mat.
     - For each goal location in  world.xGoal:
     - Uses the function Sphereworld.plot to plot the world in a first figure.
     - Sets  planner_parameters['U'] to the negative of  Total.grad.
     - it:grad-handle Calls the function Potential.planner with the problem data and the input
    arguments. The function needs to be called five times, using each one of the initial locations
    given in  x_start (also provided in !70!DarkSeaGreen2 sphereworld.mat).
     - it:plot-plan After each call, plot the resulting trajectory superimposed to the world in the
    first subplot; in a second subplot, show  u_path (using the same color and using the  semilogy
    command).
        """
        sphereWorld = SphereWorld()
        x_goal = sphereWorld.x_goal
        x_start = sphereWorld.x_start
        for goal_idx in range(x_goal.shape[1]):
            fig = plt.figure()
            plt.rcParams.update({'font.size': 8})
            ax1 = fig.add_subplot(1,2,1)
            ax1.set_title('Rep. Weight: 1, Epsilon: 0.02')
            sphereWorld.plot()
            ax2 = fig.add_subplot(1,2,2)
            ax2.set_title('Log of Potential, U')
            ax2.set_xlabel('# of Steps')
            world = sphereWorld.world
            potential = {"x_goal": x_goal[:,goal_idx].reshape(2,1),"repulsiveWeight":1,"shape":"quadratic"}
            total = Total(world,potential)
            planned_parameters = {'U': total.eval, 'control': total.grad, 'epsilon': .02, 'nb_steps':100}
            planner = Planner()
            # Run 5 times (for each x_start)
            for xstart_idx in range(x_start.shape[1]):
                x_path, u_path = planner.run(x_start[:,xstart_idx].reshape(2,1), world, potential, planned_parameters)
                ax1.scatter(x_path[0,:],x_path[1,:])
                ax2.semilogy(range(planned_parameters['nb_steps']+1),u_path)


# Planner().run_plot()
