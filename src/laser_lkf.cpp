#include "laser_lkf.h"
#include <iostream>


LaserLKF::LaserLKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in) :
  LinearKF(Mu_in, Sigma_in, t_in),
  noise_ax_(-1),
  noise_ay_(-1)
{
}


void LaserLKF::Initialize(){

  noise_ax_ = 9;
  noise_ay_ = 9;

  MatrixXd At(4,4);
  At << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;

  MatrixXd Rt(4,4);

  MatrixXd Ct(2,4);
  Ct << 1,0,0,0,
        0,1,0,0;

  MatrixXd Qt(2,2);
  Qt << 0.0225,     0,
             0,0.0225;

  LinearKF::Initialize(At, Rt, Ct, Qt);
}


void LaserLKF::Step(MeasurementPackage &meas_in){

  if (Sensor_ != meas_in.sensor_type_) return;

  double dt = meas_in.timestamp_ - previous_time_;
  previous_time_ = meas_in.timestamp_;

  dt *= 0.000001; // scale to seconds
  At_(0,2) = dt;
  At_(1,3) = dt;

  float dt_2 = dt*dt;
  float dt_3 = dt_2*dt;
  float dt_4 = dt_3*dt;
  Rt_ << dt_4/4*noise_ax_, 0, dt_3/2*noise_ax_, 0,
         0, dt_4/4*noise_ay_, 0, dt_3/2*noise_ay_,
         dt_3/2*noise_ax_, 0, dt_2*noise_ax_, 0,
         0, dt_3/2*noise_ay_, 0, dt_2*noise_ay_;

  LinearKF::Step(meas_in);
}


