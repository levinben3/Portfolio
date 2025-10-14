## High Voltage Projects
This includes the main vehicle low voltage code, which allows the vehicle to pre-charge, drive, control different sensors, and more. This vital code is written in Simulink to seamlessly integrate with any virtual testing done on code developed for the main code. This also has any other miscellaneous low voltage projects currently being developed for the vehicle. 

Any .pdf's of the same name as a .slx Simulink file are visualizations of the simulation for ease of understanding without opening in Simulink.

| File            | Explanation                                                                |
| ----------------- | ------------------------------------------------------------------ |
| regen | Initial regenerative braking control system, employing a state machine to ensure regenerative braking is only active in certain, safe scenarios. The controller uses load transfer and brake balancing to determine the optimal regenrative braking force, while preventing any dangerous effects from negative torque commands.  |
| Kalman_Filter_SOC_hv | Initial high voltage battery state of charge estimator done via a Kalman filter employing HPPC tested second-order battery model resistance and capacitance constants, as well as temperature mapping for when temperature conditions change. |
