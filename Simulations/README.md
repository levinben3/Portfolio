## Simulations Projects
All vehicle simulations work is housed in this directory. Throughout the year, we conduct simulations to help with battery sizing, inverter PID tuning, controls testing, tire modeling, and more. This is a vital and ever-expanding part of the software position at Columbia FSAE. 

Any .pdf's of the same name as a .slx Simulink file are visualizations of the simulation for ease of understanding without opening in Simulink.

| File            | Explanation                                                                |
| ----------------- | ------------------------------------------------------------------ |
| driveline_model | This simscape-based powertrain model was an initial attempt to condense the vehicle to a load and powertrain, allowing us to get estimates of straightline acceleration as well as test traction control versions in a safe environment. |
| emrax208 | This more mathematically grounded powertrain model followed up on the simscape version to help tune the Cascadia Motion Inverter PID values to improve the torque output from the driver command and also test traction control models. |
| emrax208_tc | A modified version of the mathematical powertrain model made to fit with traction control versions. |
| point_mass_sim_final | A point mass simulation designed to test the "car" on all four competition events, acceleration, skidpad, autocross, and enudurance, and return the total points received pending different parameters such as mass, energy, CLA, CDA, and more. This is based entirely in time-evolved kinematics, forces, and motor and battery stats. |
| point_mass_sim_mass_endur | This variant of the point mass simulation creates a 2D map of points for each combination of varied mass and total energy used in endurance, changing the competition power limit to do this. This proved important in 2025-26 battery sizing.|
