# shin
SHortwave INterpolation using transmissivity

These functions interpolate shortwave using the transmissivity. 

*For instantaneous irradiances, use the *INST version of the code. 

*For irradiances averaged some interval, use SHIN_v2.m. This function calculates transmissivity using the average cos(SZA) over the period in question. 

The transmissivity is then interpolated to a fine time step over which the SZA can be considered constant (5 minutes) and then averaged to the desired time step. This interpolation method conserves energy when moving to a different time step in addition to accurately handling sunrise and sunset. Linear interpolation will not conserve energy and will create inaccuracies at the shoulder periods of daylight.

These functions require the time tools found: https://github.com/klapo/time_tools
