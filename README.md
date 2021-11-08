# Robot-Configuration-Space

Functions for checking if two-link-manipulator is in collision with obstacles, plans optimal motion around obstacles.

The robot can be modeled as a set of polygons with an end-effector and the boundaries of obstacles in its environment as stars. The edge of each robot is traced in red if it collides with an obstacle, green otherwise.

![Figure](https://user-images.githubusercontent.com/61487056/140775658-35227e30-9f32-4e5d-9324-2177140ed89f.png)
![Figure2](https://user-images.githubusercontent.com/61487056/140776051-7fcd56ef-c502-4fcc-bed0-182f486b7e96.png)

When planning optimal motion with gradient descent on artifical potential functions, a quick and effective model to implement for obstacles is as spheres with a "distance of influence". Improvements upon linear potential functions (left picture below) can be made with quadratic potential functions (right picture), which tends to converge quicker. The trajectories plotted in different colors are outputted from an iterative scatter plot of the points' gradient descents, from each initial position to each position after 100 iterations. The two stars in the plot represent goals that can be inputted into the algorithm and serve as global minima of the potential field.

![Figure5](https://user-images.githubusercontent.com/61487056/140776543-55c4e1a9-290a-4e01-944d-65bf669e5f46.png)
![Figure14](https://user-images.githubusercontent.com/61487056/140776127-498d4e60-47b0-4f29-bb57-873e39683ed2.png)
