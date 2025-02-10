# Timing of motion data for the Halo Photonics lidar from VectorNav and POSMV: A Report 
## gleaned from experiments in [timing_lidar.jl](timing_lidar.jl) and [lidar_turbulence.jl](lidar_turbulence.jl)

## Leg 1 
Leg 1 has VectorNav data, presumably initially
physically aligned with the coordinate system of the lidar. 
At least its heave is aligned with the vertically-staring 
lidar beams.

VectorNav clock drifts a few seconds and then corrects, 
or precesses, with a 13-hour cycle 
(previously seemed like a 51-hour cycle), compared to the POSMV. 
The time lag offset (in discrete steps) is computed by matching the phase of the VectorNav roll to the POSMV pitch. Relative time offset is computed by counting 0.05 seconds per timestep.

<!-- figure -->
not as good:
![POSMV](./Vn-POSMV_lag.png "VectorNav-POSMV lag")
seems like a better match was achieved this time:
![VectorNav-POSMV time offset](./Vn-POSMV_lag_0601-0605.png)
##### Caption: (+) Positive lag means Vn start index is shifted forward, i.e., Vn signal lags POSMV. (Data in this shifted index is (moved backward and) aligned with the unshifted POSMV start.) (-) Negative lag means negative index of Vn lines up with start of POSMV epoch, i.e., Vn signal leads POSMV.
<!-- end figure -->

This is from syncing the ship's rotation, eliding the 18 s time offset between GPS and UTC.

The POSMV clock runs slow compared to its own internal GPS clock.
A large time offset accumulates in leg 1.
In leg 2 the clock resets to the GPS quasi-regularly every ~2 days. 
We reconstruct the POSMV time axis
to minimize the difference with the GPS (red). 

VectorNav time is 18 s ahead of POSMV, verified by plotting with DateTime x-axis. 
Amplitude modulation matches on a longer time envelope, 
and specific local offsets match the phase extremely well, with remaining time offsets mostly less than 0.5 s.
![DateTime POSMV-VN](DateTime_x_PMVpitch-VNroll_eg.png)

This agrees with 18 s GPS leapseconds.
By convention GPS time is +18 s ahead of UTC in May-June 2024.
The VectorNav time is in GPS convention and the POSMV (computer and GPS) is in UTC convention, resulting in the 18 s offset due to the GPS leapseconds.

![POSMV](./POSMV_pashr-gps.png "VectorNav-POSMV lag")


### Leg 2, part 1
There is no VectorNav data until June 4, when data starts.
There is data from the ship's POSMV starting 
around 5.18. 

Synchronization of pitch and roll
is demonstrated in `vectornav.ipynb`.
