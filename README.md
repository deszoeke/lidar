# Routines for reading and processing Doppler wind lidar 
for the 2024 ASTRAL/EKAMSAT cruise.

(c) Simon de Szoeke

## Data to read
  - Halo Photonics XR:
    `User` wind profiles and
    `Stare` vertical stares
  - VectorNav accelerometer
  - Navigation: ship's navigation message or NOAA PSL  

## Turbulence processing
Estimate vertical velocity variance and tubulent kinetic energy (TKE) dissipation rate from vertical stares.
