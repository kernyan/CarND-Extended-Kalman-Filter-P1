# Extended Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

This repository utilizes a parametric kalman filter framework to fuse radar and lidar sensor. They are implemented towards Udacity's Self-Driving Car's Term 2 Simulator.


[//]: # (Image References)

[image1]: diag.png "Parametric Kalman Filter architecture"

# Kalman Filter Architecture

The Kalman Filter architecture consists of two main building blocks, 
1. Individual filters
2. Fusion of filters

The fusion of filters is governed by a controller that adds multiple filters as needed. The fusion controller holds references to a multivariate gaussian mean and sigma, which the individual filters operate on. 

Currently supported filters

1. Linear - Based on Bayes filter algorithm
2. Extended - Based on Jacobian approximation for non-linear distributions
3. Unscented - Based on sigma-points as approximation of non-linear distributions

Additionally, the filters are implemented in the context of sensor fusion with different kinematic models. The table below shows the kinematic models available

| Kalman Filter      |  Radar	       |    Lidar           |
|:-------------------|:---------------:|:-------------------|
| Linear             |  -  	       |   Constant Vx, Vy  |
| Extended    	     | Constant Vx, Vy |	-           |
| Unscented	     | Constant turnrate & velocity | Constant turnrate & velocity |


The graph below shows the inheritance structure of the various filters
![alt text][image1]

## Example usage

```
FusionKF FKF(CONSTANT_VELOCITY); // Sensor fusion controller using constant Vx, Vy model
FKF.AddLaserLKF(); // adds Lidar Linear Kalman Filter to fusion controller
FKF.AddRadarEKF(); // adds Radar Extended Kalman Filter to fusion controller

FKF.ProcessMeasurements(Input) // Measurement inputs
```

```
FusionKF FKF(CONSTANT_TURNRATE_VELOCITY); // Sensor fusion controller using CTRV
FKF.AddLaserUKF(); // adds Lidar Unscented Kalman Filter to fusion controller
FKF.AddRadarUKF(); // adds Radar Unscented Kalman Filter to fusion controller

FKF.ProcessMeasurements(Input) // Measurement inputs
```

# Input

INPUT: values provided by the simulator to the c++ program

["sensor_measurement"] => the measurement that the simulator observed (either lidar or radar)

1. For lidar: inputs are position in x and y axis
2. For radar: inputs are rho, theta, and rho dot


# Output

OUTPUT: values provided by the c++ program to the simulator

For Constant Vx and Vy kinematic model

1. ["estimate_x"] <= kalman filter estimated position x
2. ["estimate_y"] <= kalman filter estimated position y
3. ["estimate_vx"] <= kalman filter estimated velocity along x axis
4. ["estimate_vy"] <= kalman filter estimated velocity along y axis

For Constant turnrate and velocity kinematic model

1. ["estimate_x"] <= kalman filter estimated position x
2. ["estimate_y"] <= kalman filter estimated position y
3. ["estimate_v"] <= kalman filter estimated velocity
4. ["estimate_yaw"] <= kalman filter estimated yaw
5. ["estimate_yaw_dot"] <= kalman filter estimated yaw rate

Both kinematic model produces Root-Mean-Squared-Error below

1. ["rmse_x"]
2. ["rmse_y"]
3. ["rmse_vx"]
4. ["rmse_vy"]


# Compile and running the simulator

See [Udacity's Extended Kalman Filter Project page](https://github.com/udacity/CarND-Extended-Kalman-Filter-Project) on instructions to install the Simulator. Once the necessary prerequisites are completed, the main program can be built and run by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./ExtendedKF

# Reference

Kalman filter library is maintained [here](https://github.com/kernyan/KalmanFilterController)





