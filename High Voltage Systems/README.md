## High Voltage Projects
This includes the code and controls related to the high voltage systems of the vehicle. The main high voltage code controls the accumulator or high voltage battery, and the other systems handle systems related to this, such as commanding a negative torque to allow for a return of energy through regenerative braking and monitoring the state of charge and state of health of the high voltage battery.

Any .pdf's of the same name as a .slx Simulink file are visualizations of the simulation for ease of understanding without opening in Simulink.

| File            | Explanation                                                                |
| ----------------- | ------------------------------------------------------------------ |
| regen | Initial regenerative braking control system, employing a state machine to ensure regenerative braking is only active in certain, safe scenarios. The controller uses load transfer and brake balancing to determine the optimal regenrative braking force, while preventing any dangerous effects from negative torque commands.  |
| Kalman_Filter_SOC_hv | Initial high voltage battery state of charge estimator done via a Kalman filter employing HPPC tested second-order battery model resistance and capacitance constants, as well as temperature mapping for when temperature conditions change. |
| accumulator | Main code to control the high voltage battery on the vehicle. This includes the battery monitoring system, the precharge and discharge state machine, charging code, and the overall tractive system. |
