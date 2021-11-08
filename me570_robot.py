
import numpy as np
import matplotlib.pyplot as plt
import me570_geometry as gm

class TwoLink:
    """ This class was introduced in a previous homework. """
    def __init__(self):
        add_y_reflection = lambda vertices: np.hstack(
            [vertices, np.fliplr(np.diag([1, -1]).dot(vertices))])

        vertices1 = np.array([[0, 5], [-1.11, -0.511]])
        vertices1 = add_y_reflection(vertices1)
        vertices2 = np.array([[0, 3.97, 4.17, 5.38, 5.61, 4.5],
                              [-0.47, -0.5, -0.75, -0.97, -0.5, -0.313]])
        vertices2 = add_y_reflection(vertices2)
        self.Polygons = (gm.Polygon(vertices1), gm.Polygon(vertices2))

    def polygons(self):
        """
        Returns two polygons that represent the links in a simple 2-D two-link manipulator.
        """
        return self.Polygons
    def kinematic_map(self, theta):
        """
        The function returns the coordinate of the end effector, plus the vertices of the links, all
        transformed according to  _1, _2.
        """
        # Define Rotation and Translation matrices
        W_R_B1 = gm.rot2d(theta[0])
        B1_R_B2 = gm.rot2d(theta[1])
        B1_T_B2 = np.array([[5],[0]])
        # Define B2_P = point w.r.t W (global frame)
        B2_P = np.array([[5],[0]])
        # Calculating effector global coordinate
        vertex_effector_transf = np.dot(W_R_B1,np.dot(B1_R_B2,B2_P)+B1_T_B2)
        # Get polygon vertices
        Poly1_verts = self.Polygons[0].vertices
        Poly2_verts = self.Polygons[1].vertices
        # Transform and store in new polygon object
        Poly1_transformed = gm.Polygon(np.dot(W_R_B1,self.Polygons[0].vertices))
        Poly2_transformed = gm.Polygon(np.dot(W_R_B1,np.dot(B1_R_B2,\
                                            self.Polygons[1].vertices)+B1_T_B2))
        return vertex_effector_transf, Poly1_transformed, Poly2_transformed

    def plot(self, theta, color):
        """
        This function should use TwoLink.kinematic_map from the previous question together with
        the method Polygon.plot from Homework 1 to plot the manipulator.
        """
        [vertex_effector_transf, polygon1_transf,
        polygon2_transf] = self.kinematic_map(theta)
        polygon1_transf.plot(color)
        polygon2_transf.plot(color)

    def is_collision(self, theta, points):
        """
        For each specified configuration, returns  True if  any of the links of the manipulator
        collides with  any of the points, and  False otherwise. Use the function
        Polygon.is_collision to check if each link of the manipulator is in collision.
        """
        num_thetas = len(theta[0])
        flag_theta = np.zeros(num_thetas)
        for i in range(num_thetas):
            Vertex_Effector,Transf_Poly1,Transf_Poly2 = TwoLink().kinematic_map(theta[:,i])
            if Transf_Poly1.is_collision(points).any():
                flag_theta[i] = 1    
            elif Transf_Poly2.is_collision(points).any():
                flag_theta[i] = 1    
        return flag_theta

    def plot_collision(self, theta, points):
        """
        This function should:
     - Use TwoLink.is_collision for determining if each configuration is a collision or not.
     - Use TwoLink.plot to plot the manipulator for all configurations, using a red color when the
    manipulator is in collision, and green otherwise.
     - Plot the points specified by  points as black asterisks.
        """
        flag_theta = self.is_collision(theta,points)
        for theta_index, flag in enumerate(flag_theta):
            if flag == True:
                color = 'r'
            else:
                color = 'g'
            self.plot(theta[:,theta_index], color)
        ax = plt.gca()
        ax.scatter(points[0,:],points[1,:],color='k')

