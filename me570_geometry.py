"""
 Please merge the functions and classes from this file with the same file from the previous
 homework assignment
"""

import numbers
import math
import numpy as np
from matplotlib import cm, pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

class Polygon:
    """ Class for plotting, drawing, checking visibility and collision with polygons. """
    def __init__(self, vertices):
        """
        Save the input coordinates to the internal attribute  vertices.
        """
        self.vertices = vertices
    def flip(self):
        """
        Reverse the order of the vertices (i.e., transform the polygon from filled in
        to hollow and viceversa).
        """
        self.vertices = np.fliplr(self.vertices)
    def plot(self, style):
        """
        Plot the polygon using Matplotlib.
        """

        x_pos = self.vertices[0,:]
        y_pos = self.vertices[1,:]
        x_direct, y_direct = np.diff(self.vertices, axis = 1)
        # Add in last vector
        x_direct = np.append(x_direct, x_pos[0] - x_pos[-1])
        y_direct = np.append(y_direct, y_pos[0] - y_pos[-1])
        # Quiver and Fill (if counterclockwise)
        plt.quiver(x_pos,y_pos, x_direct, y_direct, angles='xy',\
                       scale_units='xy',color=style, scale=1)
        if self.is_filled():
            plt.fill(x_pos,y_pos)
        # Enforce proportional axes
        ax = plt.gca()
        ax.axis('equal')

    def is_filled(self):
        """
        Checks the ordering of the vertices, and returns whether the polygon is filled in or not.
        
        """
        len_vec = len(self.vertices[0])
        # Insert first vertices at end to get last angle in loop
        verts = np.concatenate((self.vertices,np.reshape(self.vertices[0:,0],(2,1))),axis = 1)
        # Calculate winding number
        winding_num = 0
        for i in range(len_vec):
            winding_num += (verts[0][i] - verts[0][i+1])*(verts[1][i] + verts[1][i+1])
        if winding_num >= 0:
            flag = True #counterclockwisse
        else:
            flag = False #clockwise
        return flag
    
    def is_self_occcluded(self, idx_vertex, point):
        """
        Given the corner of a polygon, checks whether a given point is self-occluded or not by
        that polygon (i.e., if it is ``inside'' the corner's cone or not). Points on boundary
        (i.e., on one of the sides of the corner) are not considered self-occluded. Note that
        to check self-occlusion, we just need a vertex index  idx_vertex. From this, one can
        obtain the corresponding  vertex, and the  vertex_prev and  vertex_next that precede
        and follow that vertex in the polygon.
        """
        # Handle if idx_vertex is first
        if idx_vertex ==  0:
            vertex = self.vertices[0:,idx_vertex].reshape(2,1)
            vertex_prev = self.vertices[0:,-1].reshape(2,1)
            vertex_next = self.vertices[0:,idx_vertex+1].reshape(2,1)
        # Handle if idx_vertex is last
        elif idx_vertex == len(self.vertices[0]) - 1:
            vertex = self.vertices[0:,idx_vertex].reshape(2,1)
            vertex_prev = self.vertices[0:,idx_vertex-1].reshape(2,1)
            vertex_next = self.vertices[0:,0].reshape(2,1)
        else:
            vertex = self.vertices[0:,idx_vertex].reshape(2,1)
            vertex_prev = self.vertices[0:,idx_vertex-1].reshape(2,1)
            vertex_next = self.vertices[0:,idx_vertex+1].reshape(2,1)
        
        # Get angle betwen edges linked to idx_vertex
        theta1 = angle(vertex,vertex_prev,vertex_next)
        # Get angle between vertex_prev and point
        theta2 = angle(vertex,vertex_prev,point)
        
        if self.is_filled() != bool(np.sign(theta1)+1):
            if np.sign(theta1) == np.sign(theta2):
                if abs(theta2)>=abs(theta1):
                    flag_point = False
                else:
                    flag_point = True
            else:
                if 2*np.pi-abs(theta1) > abs(theta2):
                    flag_point = False
                else:
                    flag_point = True
        else:
            if np.sign(theta1) == np.sign(theta2):
                if abs(theta2)>=abs(theta1):
                    flag_point = True
                else:
                    flag_point = False
            else:
                if 2*np.pi-abs(theta1) > abs(theta2):
                    flag_point = True
                else:
                    flag_point = False
        return flag_point

    def is_visible(self, idx_vertex, test_points):
        """
        Checks whether a point p is visible from a vertex v of a polygon. In order to be visible,
        two conditions need to be satisfied: enumerate  point p should not be self-occluded with
        respect to the vertex v (see Polygon.is_self_occluded). The segment p--v should not collide
        with any of the edges of the polygon (see Edge.is_collision).
        """

        nb_points = len(test_points[0])
        nb_edges = len(self.vertices[0])
        flag_points = np.zeros(nb_points)
        for i in range(nb_points):
            # Check occlusion
            if self.is_self_occcluded(idx_vertex,test_points[0:,i].reshape(2,1)):
                flag_points[i] = False
            else:
                flag_points[i] = True
                # Check collision
                vertex_x1,vertex_y1 = self.vertices[0,idx_vertex],self.vertices[1,idx_vertex]
                vertex_x2,vertex_y2 = test_points[0,i],test_points[1,i]
                Edge1= Edge(np.array([[vertex_x1,vertex_x2],[vertex_y1,vertex_y2]]))
                for j in range(nb_edges):
                    vertex_x3,vertex_y3 = self.vertices[0,j],self.vertices[1,j]
                    if j == nb_edges-1:
                        vertex_x4,vertex_y4 = self.vertices[0,0],self.vertices[1,0]
                    else:
                        vertex_x4,vertex_y4 = self.vertices[0,j+1],self.vertices[1,j+1]
                    Edge2 = Edge(np.array([[vertex_x3,vertex_x4],[vertex_y3,vertex_y4]]))
                    if Edge1.is_collision(Edge2):
                        flag_points[i] = False
        return flag_points

    def is_collision(self, test_points):
        """
        Checks whether the a point is in collsion with a polygon (that is, inside for a filled in
        polygon, and outside for a hollow polygon). In the context of this homework, this function
        is best implemented using Polygon.is_visible.
        """
        numb_points = len(test_points[0])
        numb_vertices = len(self.vertices[0])
        vis = np.zeros((numb_vertices,numb_points))
        for i in range(numb_vertices):
            vis[i,:] = self.is_visible(i,test_points)
        flag_points = ~np.array(sum(vis), dtype=bool)
        if self.is_filled():
            return flag_points
        else:
            return ~flag_points


class Grid:
    """
    A function to store the coordinates of points on a 2-D grid and evaluate arbitrary
    functions on those points.

    [This class is the same as the one from HW2]
    """
    def __init__(self, xx_grid, yy_grid):
        """
        Stores the input arguments in attributes.
        """
        self.xx_grid = xx_grid
        self.yy_grid = yy_grid

    def eval(self, fun):
        """
        This function evaluates the function  fun (which should be a function)
        on each point defined by the grid.
        """

        dim_domain = [numel(self.xx_grid), numel(self.yy_grid)]
        dim_range = [numel(fun(np.array([[0], [0]])))]
        fun_eval = np.nan * np.ones(dim_domain + dim_range)
        for idx_x in range(0, dim_domain[0]):
            for idx_y in range(0, dim_domain[1]):
                x_eval = np.array([[self.xx_grid[idx_x]],
                                   [self.yy_grid[idx_y]]])
                fun_eval[idx_x, idx_y, :] = np.reshape(fun(x_eval),
                                                       [1, 1, dim_range[0]])

        # If the last dimension is a singleton, remove it
        if dim_range == [1]:
            fun_eval = np.reshape(fun_eval, dim_domain)

        return fun_eval

    def mesh(self):
        """
        Shorhand for calling meshgrid on the points of the grid
        """

        return np.meshgrid(self.xx_grid, self.yy_grid)



class Torus:
    """
    A class that holds functions to compute the embedding and display a torus and curves on it.
    """
    def phi(self, theta):
        """
        Implements equation (eq:chartTorus).
        """
        #### Input is 2*nbpoints. Output is 3*nbpoints.
        nb_points = theta.shape[1]
        x_torus = np.zeros((3,nb_points))
        for i in range(nb_points):
            [theta1, theta2] = theta[:,i].flatten()
            diag_theta2 = np.array([[np.cos(theta2),-np.sin(theta2),0],\
                                        [np.sin(theta2),np.cos(theta2),0],\
                                            [0,0,1]]
                                       )
            map_xz_plane = np.array([[1,0],[0,0],[0,1]])
            radius_translation = np.array([[3],[0],[0]])
            R_col = np.array([[np.cos(theta1)],[np.sin(theta1)]])
            x_torus[:,i] = np.dot(diag_theta2,np.dot(map_xz_plane,R_col)+radius_translation).flatten()
        return x_torus

    def plot_charts(self):
        """
        For each one of the chart domains U_i from the previous question:
        - Fill a  grid structure with fields  xx_grid and  yy_grid that define a grid of regular
          point in U_i. Use nb_grid=33.
        - Call the function Grid.eval with argument Torus.phi.
        - Plots the surface described by the previous step using the the Matplotlib function
        ax.plot_surface (where  ax represents the axes of the current figure) in a separate figure.
        Plot a final additional figure showing all the charts at the same time.   To better show
        the overlap between the charts, you can use different colors each one of them,
        and making them slightly transparent.
        """
        nb_grid = 33
        # Plot first chart.
        xx_U1 = np.linspace(0-.1,np.pi+.1,num=nb_grid)
        yy_U1 = np.linspace(0-.1,2*np.pi+.1,num=nb_grid)
        grid_U1 = Grid(xx_U1,yy_U1)
        grid_U1_eval = grid_U1.eval(self.phi)
        X1 = grid_U1_eval[:,:,0]
        Y1 = grid_U1_eval[:,:,1]
        Z1 = grid_U1_eval[:,:,2]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X1, Y1, Z1,color = 'b',alpha=1)
        #Plot second chart
        xx_U2 = np.linspace(np.pi+.1,2*np.pi+.1,num=nb_grid)
        yy_U2 = np.linspace(0-.1,2*np.pi+.1,num=nb_grid)
        grid_U2 = Grid(xx_U2,yy_U2)
        grid_U2_eval = grid_U2.eval(self.phi)
        X2 = grid_U2_eval[:,:,0]
        Y2 = grid_U2_eval[:,:,1]
        Z2 = grid_U2_eval[:,:,2]
        ax.plot_surface(X2, Y2, Z2,color = 'r',alpha=.75)
        
        plt.show()


    def phi_push_curve(self, a_line, b_line):
        """
        This function evaluates the curve x(t)= phi_torus ( phi(t) )  R^3 at  nb_points=31 points
        generated along the curve phi(t) using line_linspaceLine.linspace with  tMin=0 and  tMax=1,
        and a, b as given in the input arguments.
        """
        t_min=0
        t_max=1
        nb_points = 20
        x_points = np.zeros((3,nb_points))
        theta = line_linspace(a_line, b_line, t_min, t_max, nb_points)
        x_points = self.phi(theta)
        return x_points

    def plot_charts_curves(self):
        """
        The function should iterate over the following four curves:
        - 3/4*pi0
        - 3/4*pi3/4*pi
        - -3/4*pi3/4*pi
        - 03/4*pi  and  b=np.array([[-1],[-1]]).
        The function should show an overlay containing:
        - The output of Torus.plotCharts;
        - The output of the functions torus_pushCurveTorus.pushCurve for each one of the curves.
        """
        a_line = np.array([[3/4*np.pi],[0]])
        b_line = np.array([[-1],[0]])
        x_points = self.phi_push_curve(a_line,b_line)
        ax = plt.gca(projection='3d')
        ax.plot(x_points[:,0], x_points[:,1], x_points[:,2], 'green')
        self.plot_charts()
    
class Edge:
    """ Class for storing edges and checking collisions among them. """
    def __init__(self, vertices):
        """
        Save the input coordinates to the internal attribute  vertices.
        """
        self.vertices = vertices

    def is_collision(self, edge):
        """
        Returns  True if the two edges intersect.  Note: if the two edges overlap but are colinear,
        or they overlap only at a single endpoint, they are not considered as intersecting (i.e.,
        in these cases the function returns  False). If one of the two edges has zero length, the
        function should always return the result that edges are non-intersecting.
        """
        # tolerance
        t = -0.01
        # Break up edges into points
        x1,y1 = self.vertices[0,0],self.vertices[1,0]
        x2,y2 = self.vertices[0,1],self.vertices[1,1]
        x3,y3 = edge.vertices[0,0],edge.vertices[1,0]
        x4,y4 = edge.vertices[0,1],edge.vertices[1,1]
        # Get slope y-intercept, and domains
        if x4==x3:
            m2 = np.power(10,10)
        else:
            m2 = (y4-y3)/(x4-x3)
        if x2==x1:
            m1 = np.power(10,10)
        else:
            m1 = (y2-y1)/(x2-x1)
        # False if colinear
        if m1 == m2:
            return False 
        # False if one of the endpoints equals another
        elif (x1 == x3 and y1 == y3) or (x1 == x4 and y1 == y4)\
            or (x2  == x3 and y2 == y3) \
            or (x2 == x4 and y2 == y4):
            return False
        # False if length of an edge is zero
        elif np.linalg.norm(self.vertices) == 0 or np.linalg.norm(edge.vertices) == 0:
            return False
        else:
            # Get y-intercept for each edge
            b1 = y1 - m1*x1
            b2 = y3 - m2*x3
            # Get solution by setting equations equal
            x_solution = (b2-b1)/(m1-m2)
            y_solution = m1*x_solution + b1
            # Check if solution is in domain and range of edges
            domain1,domain2 = np.array([x1,x2]),np.array([x3,x4])
            range1,range2 = np.array([y1,y2]),np.array([y3,y4])
            domain1.sort(); domain2.sort(); range1.sort(); range2.sort()
            # Create bool from domains
            bool_domain_intersection = [x_solution-domain1[0]>=t, domain1[1]-x_solution>=t,\
                                        x_solution-domain2[0]>=t, domain2[1]-x_solution>=t]
            # Create bool from ranges
            bool_range_intersection = [y_solution-range1[0]>=t, range1[1]-y_solution>=t,\
                                       y_solution-range2[0]>=t, range2[1]-y_solution>=t]
            if all(bool_range_intersection + bool_domain_intersection):
                return True
        return False
    
class Sphere:
    """ Class for plotting and computing distances to spheres (circles, in 2-D). """
    def __init__(self, center, radius, distance_influence):
        """
        Save the parameters describing the sphere as internal attributes.
        """
        self.center = center
        self.radius = radius
        self.distance_influence = distance_influence

    def plot(self, color):
        """
        This function draws the sphere (i.e., a circle) of the given radius, and the specified color,
    and then draws another circle in gray with radius equal to the distance of influence.
        """
        # Get current axes
        ax = plt.gca()
        # Add circle as a patch
        if self.radius > 0:
            # Circle is filled in
            kwargs = {'facecolor': (0.3, 0.3, 0.3)}
            radius_influence = self.radius + self.distance_influence
        else:
            # Circle is hollow
            kwargs = {'fill': False}
            radius_influence = -self.radius - self.distance_influence

        center = (self.center[0, 0], self.center[1, 0])
        ax.add_patch(
            plt.Circle(center,
                       radius=abs(self.radius),
                       edgecolor=color,
                       **kwargs))

        ax.add_patch(
            plt.Circle(center,
                       radius=radius_influence,
                       edgecolor=(0.7, 0.7, 0.7),
                       fill=False))
        ax.axis("equal")

    def distance(self, points):
        """
        Computes the signed distance between points and the sphere, while taking into account whether
    the sphere is hollow or filled in.
        """
        sign_rad = np.sign(self.radius)
        if len(points.shape) == 1:
            points = points.reshape(2,1)
        d_points_sphere = np.array([])
        for pnt_indx in range(points.shape[1]):
            point = points[:,pnt_indx].reshape(2,1)
            pnt_dist = np.linalg.norm(point - self.center)-abs(self.radius)
            d_points_sphere = np.append(d_points_sphere,pnt_dist*sign_rad)
        
        return d_points_sphere

    def distance_grad(self, sphere, points):
        """
        Computes the gradient of the signed distance between points and the sphere, consistently with
    the definition of Sphere.distance.
        """
        # print('POINTS')
        # print(points)
        # print('whatttt')
        if len(points.shape) == 1:
            points = points.reshape(2,1)
        
        grad_d_points_sphere = np.array([[],[]])
        for pnt_indx in range(points.shape[1]):
            point = points[:,pnt_indx].reshape(2,1)
            if (point == sphere.center).all():
                dist_grad = np.zeroz((2,1))
            else:
                ####
                # vec_pt_center = point - sphere.center
                # theta = np.arctan(vec_pt_center[1,0]/vec_pt_center[0,0])
                # vec_dist = sphere.distance(point)[0]
                # x_vec_dist = vec_dist * np.cos(theta)
                # y_vec_dist = vec_dist * np.sin(theta)
                # vec_sphere_point = np.array([[x_vec_dist],[y_vec_dist]])
                # dist_grad = vec_sphere_point/vec_dist
                # dist = self.distance(points)
                ####
                dist = self.distance(point)[0]
                vec_pt_center = (point - sphere.center)/np.linalg.norm(point - sphere.center)
                vec = sphere.center + sphere.radius*vec_pt_center
                dist_grad = (point-vec)/dist
            grad_d_points_sphere = np.append(grad_d_points_sphere,dist_grad,axis=1)
        # print(grad_d_points_sphere)
        return grad_d_points_sphere
    
def angle(vertex0, vertex1, vertex2, angle_type='signed'):
    """
    Compute the angle between two edges  vertex0-vertex1 and  vertex0-vertex2 having an endpoint in
    common. The angle is computed by starting from the edge  vertex0-- vertex1, and then
    ``walking'' in a counterclockwise manner until the edge  vertex0-vertex2 is found.
    The angle is computed by starting from the vertex0-vertex1 edge, and then “walking” in a
    counterclockwise manner until the is found.
    """
    # tolerance to check for coincident points
    tol = 2.22e-16
    # compute vectors corresponding to the two edges, and normalize
    vec1 = vertex1 - vertex0
    vec2 = vertex2 - vertex0
    norm_vec1 = np.linalg.norm(vec1)
    norm_vec2 = np.linalg.norm(vec2)
    if norm_vec1 < tol or norm_vec2 < tol:
        # vertex1 or vertex2 coincides with vertex0, abort
        edge_angle = math.nan
        return edge_angle
    vec1 = vec1 / norm_vec1
    vec2 = vec2 / norm_vec2
    # Transform vec1 and vec2 into flat 3-D vectors,
    # so that they can be used with np.inner and np.cross
    vec1flat = np.vstack([vec1, 0]).flatten()
    vec2flat = np.vstack([vec2, 0]).flatten()
    c_angle = np.inner(vec1flat, vec2flat)
    if c_angle == 0:
        c_angle = np.pi
    s_angle = np.inner(np.array([0, 0, 1]), np.cross(vec1flat, vec2flat))
    edge_angle = math.atan2(s_angle, c_angle)
    angle_type = angle_type.lower()
    if angle_type == 'signed':
        # nothing to do
        pass
    elif angle_type == 'unsigned':
        edge_angle = math.modf(edge_angle + 2 * math.pi, 2 * math.pi)
    else:
        raise ValueError('Invalid argument angle_type')
    return edge_angle



def numel(var):
    """
    Counts the number of entries in a numpy array, or returns 1 for fundamental numerical
    types
    """
    if isinstance(var, numbers.Number):
        size = int(1)
    elif isinstance(var, np.ndarray):
        size = var.size
    else:
        raise NotImplementedError(f'number of elements for type {type(var)}')
    return size


def rot2d(theta):
    """
    Create a 2-D rotation matrix from the angle  theta according to (1).
    """
    rot_theta = np.array([[math.cos(theta), -math.sin(theta)],
                          [math.sin(theta), math.cos(theta)]])
    return rot_theta


def line_linspace(a_line, b_line, t_min, t_max, nb_points):
    """
    Generates a discrete number of  nb_points points along the curve
    (t)=( a(1)  t+ b(1), a(2) t+b(2))  R^2 for t ranging from  tMin to  tMax.
    """
    beginning_pt = t_min*a_line + b_line
    end_pt = t_max*a_line + b_line
    theta_points = np.linspace(beginning_pt,end_pt,nb_points)
    return theta_points[:,:,0].T


def clip(val, threshold):
    """
    If val is a scalar, threshold its value; if it is a vector, normalized it
    """
    if isinstance(val, np.ndarray):
        val_norm = np.linalg.norm(val)
        if val_norm > threshold:
            val /= val_norm
    elif isinstance(val, numbers.Number):
        if val > threshold:
            val = threshold
    else:
        raise ValueError('Numeric format not recognized')
    return val


def field_plot_threshold(f_handle, threshold=10, nb_grid=61):
    """
    The function evaluates the function  f_handle on points placed on a regular grid.
    """

    xx_grid = np.linspace(-11, 11, nb_grid)
    yy_grid = np.linspace(-11, 11, nb_grid)
    grid = Grid(xx_grid, yy_grid)

    f_handle_clip = lambda val: clip(f_handle(val), threshold)
    f_eval = grid.eval(f_handle_clip)

    [xx_mesh, yy_mesh] = grid.mesh()
    f_dim = numel(f_handle_clip(np.zeros((2, 1))))
    if f_dim == 1:
        # scalar field

        fig = plt.gcf()
        
        axis = fig.add_subplot(111, projection='3d')

        axis.plot_surface(xx_mesh,
                          yy_mesh,
                          f_eval.transpose(),
                          cmap=cm.gnuplot2)
        axis.view_init(90, -90)
    elif f_dim == 2:
        # vector field
        # grid.eval gives the result transposed with respect to what meshgrid expects
        f_eval = f_eval.transpose((1, 0, 2))
        # vector field
        plt.quiver(xx_mesh,
                   yy_mesh,
                   f_eval[:, :, 0],
                   f_eval[:, :, 1],
                   angles='xy',
                   scale_units='xy')
    else:
        raise NotImplementedError(
            'Field plotting for dimension greater than two not implemented')
    
    plt.xlabel('x')
    plt.ylabel('y')


# center = np.array([[3.75],[0]])
# radius = 1
# distance_influence = 10
# Sphre = Sphere(center, radius, distance_influence)
# # pt = np.array([[3.76072202],[3.22443907]])
# pt = np.array([[3.7499],[3.22443055]])
# print('dist')
# print(Sphre.distance(pt))
# print('grad_dist')
# print(Sphre.distance_grad(Sphre,pt))

# sphere=Sphere(2*np.ones((2,1)),-2,3)
# f_handle= lambda point: sphere.distance(point)
# field_plot_threshold(f_handle,10)
