# eclipse_finding
This project aims to provide general purpose code to find eclipses using Skyfield.

# Background

Eclipses involve three bodies, the eclipsee (A), the eclipser (B), and the screen (C). The light from the eclipsee is (partially) blocked by the eclipser before it reaches the screen. In case of a solar eclipse, the Sun is the eclipsee, the Moon is the eclipser, and the Earth is the screen onto which the shadow of the Moon falls. During a lunar eclipse the light from the Sun (eclipsee) is blocked by the Earth (eclipser) before it falls on the Moon (screen). Lastly, the eclipse is observed from somewhere (observer) - this is usually the Earth (both in case of the solar and lunar eclipse) but could potentially by a satellite, or any other body.
![coordinate_system](https://github.com/user-attachments/assets/76f8e22f-5745-4059-be8f-5a0602e98e6d)

Typically, the eclipsee is a star, such that the shadow of the eclipser falls onto the screen. There are three types of shadow: the umbra, the penumbra, and the antumbra. Within the umbra the eclipser completely covers the light source (eclipsee). The antumbra extends beyond the tip of the umbra, within which the eclipser is completely in front of the light source (eclipsee) but too small to completely cover it. The penumbra describes the part of the shadow within which the eclipser is only partially in front of the light source (eclipsee).

In order to be able to quickly calculate whether an eclipse occurs at time $t$, we need the positions of the three bodies A, B, and C in a common coordinate system - in fact, we need the apparent position of A and B, seen from C. Let's call these position vectors $\vec{a}$, $\vec{b}$, and $\vec{c}$ - see Figure above. The line connecting the (bary)centers of A and B is $\vec{d} = \vec{b} - \vec{a}$. Checking whether the line $\vec{d}$, extended beyond $B$, intersects a sphere with the radius of $C$ ($r_C$) centered on $\vec{c}$, provides a quick way to check whether an eclipse occurred. The line connecting $A$ and $B$ is given by

$$\vec{l}(s) = \vec{a} + s\left ( \vec{b} - \vec{a} \right ),$$

and to check for a intersection with a sphere of radius $r_C$ around $\vec{c}$, we need to solve for the parameter $s$ such that

$$\left ( \vec{l}(s) - \vec{c} \right ) \cdot \left ( \vec{l}(s) - \vec{c} \right ) = r_C^2, $$

or

$$\left ( \vec{a} + s\left ( \vec{b} - \vec{a} \right ) - \vec{c} \right ) \cdot \left ( \vec{a} + s\left ( \vec{b} - \vec{a} \right ) - \vec{c} \right ) = r_C^2.$$

Rearrainging we get

$$s^2 \left ( \vec{b} - \vec{a} \right ) \cdot \left ( \vec{b} - \vec{a} \right ) + 2s \left ( \vec{b} - \vec{a} \right ) \cdot \left ( \vec{a} - \vec{c} \right ) + \left ( \vec{a} - \vec{c} \right ) \cdot \left ( \vec{a} - \vec{c} \right ) - r_C^2 = 0$$

which is a quadratic equation of the form $as^2 + bs + c = 0$,

$$a = \left ( \vec{b} - \vec{a} \right ) \cdot \left ( \vec{b} - \vec{a} \right )$$

$$b = 2 \left ( \vec{b} - \vec{a} \right ) \cdot \left ( \vec{a} - \vec{c} \right )$$

$$c = \left ( \vec{a} - \vec{c} \right ) \cdot \left ( \vec{a} - \vec{c} \right ) - r_C^2$$

and the two solutions

$$s_{1,2} = \dfrac{-b \pm \sqrt{b^2 - 4ac}}{2a}.$$

Only if $b^2 \ge 4ac$ is there a solution. As an aside, we could assume body $C$ to be an ellipsoid, which would make the calculations slightly more involved - for simplicity sake, we here assume all bodies are spheres. Once we find the solutions $s_{1,2}$, we can find the intersection points, and select the one closer to $B$. If we wanted to express this location in some local coordinate system (like latitude/longitude), we would have to transform the location of the intersection accordingly.

Of course, the line connecting the centers does not account for the spatial extent of the umbra, and we could miss cases in which the centerline $\vec{l}(s)$ misses the screen $C$, but some part of the extended umbra (and penumbra) still cause an eclipse on $C$.

# Spatial extent of eclipses

In order to estimate the extent of the shadow, we can use the following raytracing strategy. We can find a vector that is perpendicular to $\vec{d}$ (and $\vec{a}$) by calculating $\vec{r} = \vec{a} \times \vec{d}$ (see Figure above). If we scale the length of $\vec{r}$ to the radius of $A$, $r_A$, we have a radius vector for $A$ that, if rotated around $\vec{d}$ will trace a line on the rim of $A$, as seen from $B$. We can use the same vector $\vec{r}$, but centered at $\vec{b}$ and scaled to $r_B$ to trace along the rim of $B$ when rotated around $\vec{d}$. For each point on the rim of $A$, we trace lines (or rays) to all the points on the rim of $B$, and check whether the extended line crosses $C$ (a sphere of radius $r_C$ centered at $\vec{c}$). From the angle the line connecting two points on the rim makes with the centerline, we can find out whether that particular line is in the umbra or penumbra.

The line from one point on the rim of $A$ to another point on the rim of $B$ is given by

$$\vec{l}(s) = \vec{a} + R_A \vec{r} + s \left (\vec{b} + R_B\vec{r} - \vec{a} - R_A\vec{r} \right ),$$

where $R_A$ and $R_B$ are matrices facilitating the scaling of $\vec{r}$, and its rotation around $\vec{d}$. Assuming that the length of the radius vector is $r^\prime$ and the body's radius is $r$, the rotation matrix is given by [^1]

$$R(\vec{d}, \theta) = \begin{pmatrix}
		k u_x^2 (1 - \cos{\theta}) + \cos{\theta} & u_x u_y (1 - \cos{\theta}) - u_z \sin{\theta} & u_x u_z (1 - \cos{\theta}) + u_y \sin{\theta} \\
		(u_x u_y (1 - \cos{\theta}) + u_z \sin{\theta} & k u_y^2 (1 - \cos{\theta}) + \cos{\theta} & u_y u_z (1 - \cos{\theta}) - u_x \sin{\theta} \\
		(u_x u_z (1 - \cos{\theta}) - u_y \sin{\theta} & u_y u_z (1 - \cos{\theta}) + u_x \sin{\theta} & k u_z^2 (1 - \cos{\theta}) + \cos{\theta})
		\end{pmatrix},$$

where $\theta$ is some (arbitrary) angle between 0 and 2 $\pi$ and $\vec{u} = \vec{d}/d$ and $k = r/r^\prime$. For each combination of $R_A$ and $R_B$ we calculate (if they exist) $s_{1,2}$, analogously to before. The angle $a$ between the rimlines and the centerline is

$$\cos{a} = \dfrac{\left ( \vec{b} - \vec{a} \right ) \cdot \left ( \vec{b} + R_B\vec{r} - \vec{a} - R_A\vec{r} \right )}{\left | \vec{b} - \vec{a} \right | \left | \vec{b} + R_B\vec{r} - \vec{a} - R_A\vec{r} \right |}.$$

![eclipse](https://github.com/user-attachments/assets/6375287b-f6d1-4530-9747-51a9722edd60)

In the Figure above, different angles for the umbra, penumbra, and antumbra are shown. Here, angles are positive clockwise - in the second sketch, $a_1$ is therefore negative. When the angle $a$ between the rimline and the centerline is greater than $a_1$, but less than $a_3$, the ray falls into the umbra. If $a$ is greater than both $a_1$ and $a_3$, but smaller than $a_2$, it falls into the penumbra. The angles are defined as

$$\tan{a_1} = \dfrac{r_A - r_B}{d}$$

$$\tan{a_2} = \dfrac{r_A + r_B}{d}$$

$$\tan{a_3} = \dfrac{r_A + r_B - d_2 \tan{a_1}}{d + d_2}$$

If $d_2$ is longer than $d_3 = \dfrac{r_B}{\tan{a_1}}$, we are in the antumbra case, which is weird.

Once we have the coordinates of the intercepts between the rimlines and the screen $C$, we can transform these into a local coordinate system like a planetrary latitude/longitude system.

[^1]: Rotation matrix on Wikipedia(https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle)
