## Low Voltage Projects
This includes the main vehicle low voltage code, which allows the vehicle to pre-charge, drive, control different sensors, and more. This vital code is written in Simulink to seamlessly integrate with any virtual testing done on code developed for the main code. This also has any other miscellaneous low voltage projects currently being developed for the vehicle. 

Any .pdf's of the same name as a .slx Simulink file are visualizations of the simulation for ease of understanding without opening in Simulink.

| File            | Explanation                                                                |
| ----------------- | ------------------------------------------------------------------ |
| main | Main low voltage vehicle code which includes the vehcile state machine, sensors, pedals curve, traction control, vehicle dash code, and other vital safety systems. |
| Kalman_Filter_SOC | Initial low voltage battery state of charge estimator done via a Kalman filter. |
