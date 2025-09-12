## Traction Control Project
This 2024-25 project focused on developing a forward speed kalman filter to help calculate slip ratio and then a feed-forward traction control system. This traction control system was initially employed for the straightline acceleration event for the Formula SAE competition, but we are now implementing it to run on curves as well. 

Any .pdf's of the same name as a .slx Simulink file are visualizations of the simulation for ease of understanding without opening in Simulink.

| File            | Explanation                                                                |
| ----------------- | ------------------------------------------------------------------ |
| traction_control_testing | This file contains data feeding into the kalman filter for forward speed and multiple different traction controllers. The one currently used is the same as the one in the main LV code, but there are other simpler versions used at different points, such as a simpler, first order model for the endurance event as an initial start. The controller has an open loop and closed loop controller design with two torques calculated and added together. |
| Kalman_Filter_Mod | This is the object-oriented Matlab system that runs the Kalman filter for forward speed using wheel speeds and the accelerometer. |
