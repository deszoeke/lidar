# Compute TKE dissipation rate from lidar quasi vertical stares

The stabilization table was rolling with the ship, so it sampled some component of the mean horizontal wind.

## The ship coordinate system
Define the ship coordinate to be "northward" toward the bow and "eastward" to starboard. The direction a wind is from is 0 from north, and increases to 90° easterly. The $x$ coordinate points forward (north) in this coordinate system, and the $y$ coordinate points to starboard (east). $z$ is upward for a right-handed coordinate system.

Pitch $\theta$ (y-axis) and roll $\phi$ (x axis) rotations are defined in this coordinate system. Positive roll angle is a rotation of the port upward; positive pitch angle is a rotation of the bow upward.

## Winds
The vertical wind is $w$. The ship-relative wind 
is $W$.

Ship-relative horizontal wind is
$ (U, V) = $ (-speed cos(dir), -speed sin(dir)).

Radial velocity away from the lidar is $\hat{w}$ with
$$
\hat{w}^2 = (-U sin\theta)^2 
+ (V cos\theta sin\phi)^2 
+ (W cos\theta cos\phi)^2 
$$.

## Relative velocity
We solve the measured Doppler radial velocity $\hat w$ for the air velocity,
$$
W^2 = \frac{ \hat{w}^2 
- (U sin\theta)^2 
- (V cos\theta sin\phi)^2 }
{(cos\theta cos\phi)^2}
$$.

## Vectors among the lidar, and two target volumes
The radial vector to the target volume is
$$
{\bf \hat{r}} =
\begin{pmatrix}
\hat x \\ \hat y \\ \hat z
\end{pmatrix}
=
({\rm range})
\begin{pmatrix}
-sin\theta \\
 cos\theta sin\phi \\
 cos\phi cos\theta
\end{pmatrix}
$$

The displacement between two *simultaneous* target
volumes is
$$
\bf{\hat r_1 - \hat r_2} = 
\begin{pmatrix}
\hat x \\
\hat y \\
\hat z
\end{pmatrix}
=
{\rm range}
\begin{pmatrix}
-sin \theta_1-(-sin \theta_2) \\
cos\theta_2 sin \phi_1 - cos\theta_2 sin \phi_2 \\
cos\theta_1 cos \phi_1 - cos\theta_2 cos\phi_2
\end{pmatrix}
$$

In the time $\tau$ the ship moves relative to the air a horizontal displacement of $(X, Y) = \tau(U,V)$. The total displacement of the sample volumes is 
$$
{\bf r_2 - r_1} = 
\begin{pmatrix}
\hat x + X\\
\hat y + Y\\
\hat z
\end{pmatrix}
$$
Expanding, the distance between the sample volumes is 
$$
R = |{\bf r_2 - r_1}| =
\begin{pmatrix}{|}
\tau U - &{\rm range}(sin\theta_1 + sin\theta_2) \\
\tau V + &{\rm range}(cos\theta_1 sin\phi_1 - cos\theta_2 sin\phi_2) \\
 & {\rm range}(cos\theta_1 cos\phi_1 - cos\theta_2 cos\phi_2)
\end{pmatrix}
$$
This is the displacement argument of the (second order) structure function 
$$
D(R) = \overline{(w_2-w_1)^2}
$$,
which behaves as $D = N + AR^{2/3}$.
The structure function is averaged as a function of $R$. We will retain sample pairs of $(D,R)$ and bin average in $R$, perhaps with bins equally populated by sample pairs.

## Ship motion fluctuations
The ship heave, pitch, and roll also induce fluctuating velocities in the transmitter that must be added. (The separation of scale between the mean ship-relative velocity and these fast fluctuations breaks down in maneuvers where the ship changes direction quickly.) Heave, pitch, and roll is either measured at the lidar on the stabilization table, where it can be added directly to the radial Doppler $\hat w$, measured in the ship orientation at the location of the lidar, or at a reference location on the ship (see survey).

Presently, I read a Notre Dame VectorNav file for measuring and correcting the motion of the table, but I don't know the quantities and units. The ship's POSMV IMU is also available from a reference location. If we use the POSMV, we need to calculate the velocities and moments from the heave and angular rates.

Assuming the ship is 
rigid, the heave is the same everywhere. Angular rates induce motion in proportion to the moment
${\bf L} = (L_x, L_y, L_z)$ 
be the displacement vector from the reference location 
(of the POSMV) to the lidar. At the lidar the ship
moves by 
$$
\dot{x} = -\dot\theta L_z\\
\dot{y} = \dot\phi L_z\\
\dot{z} = \dot\theta L_x - \dot\phi L_y
$$.
Varying yaw is ignored.

