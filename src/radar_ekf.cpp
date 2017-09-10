#include "radar_ekf.h"
#include "tools.h"
#include <iostream>


RadarEKF::RadarEKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in) :
  ExtendedKF(Mu_in, Sigma_in, t_in),
  noise_ax_ (-1),
  noise_ay_ (-1)
{
}


void RadarEKF::Step(MeasurementPackage &meas_in){

  if (Sensor_ != meas_in.sensor_type_) return;

  double dt = meas_in.timestamp_ - previous_time_;
  previous_time_ = meas_in.timestamp_;

  dt *= 0.000001; // scale to seconds
  G_(0,2) = dt;
  G_(1,3) = dt;

  auto G_in = G_;
  g_ = [G_in](VectorXd Mu_in){ // prediction function is assigned here because depends on dt
    return G_in*Mu_in;
  };

  float dt_2 = dt*dt;
  float dt_3 = dt_2*dt;
  float dt_4 = dt_3*dt;
  Rt_ << dt_4/4*noise_ax_, 0, dt_3/2*noise_ax_, 0,
         0, dt_4/4*noise_ay_, 0, dt_3/2*noise_ay_,
         dt_3/2*noise_ax_, 0, dt_2*noise_ax_, 0,
         0, dt_3/2*noise_ay_, 0, dt_2*noise_ay_;

  CalculateMuBar();
  CalculateSigmaBar();
  MatrixXd S = CalculateMeasurementVar();
  MatrixXd Ht = H_(MuBar_);
  MatrixXd K = SigmaBar_ * Ht.transpose() * S.inverse();
  VectorXd zhat = CalculatePredictedMeasurement();
  VectorXd y = meas_in.raw_measurements_ - zhat;

  NormToPi(1, y); // normalize theta to between +/- pi

  if (meas_in.raw_measurements_(1) * zhat(1) <= 0){
    y << 0,0,0; // do not update when theta flips sign
  }
  Mu_ = MuBar_ + K * y;
  Sigma_ = SigmaBar_ - K * S * K.transpose();
}


void RadarEKF::Initialize(){

  noise_ax_ = 9;
  noise_ay_ = 9;

  MatrixXd Gt(4,4);
  Gt << 1,0,0,0,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1;

  function<VectorXd (VectorXd)> g_in;
  MatrixXd Rt(4,4);

  function<VectorXd (VectorXd)> h_in = [](VectorXd MuBar_in){
    VectorXd zhat(3);

    float px = MuBar_in(0);
    float py = MuBar_in(1);
    float vx = MuBar_in(2);
    float vy = MuBar_in(3);

    float p = sqrt(px*px + py*py);
    float phi = atan2(py, px);
    if (fabs(py) < 0.0001 && fabs(px) < 0.0001) cout << "Warning: atan2(0,0) in RadarEKF::Initialize encountered\n";

    float p_dot = 0;
    if (fabs(p) < 0.0001) {
      p = 0.0001; // arbitrary small number
      cout << "Warning: div0 in RadarEKF::Initialize encountered\n";
    }

    p_dot = (px*vx+py*vy)/p;

    zhat << p, phi, p_dot;
    return zhat;
  };

  function<MatrixXd (VectorXd)> H_in = [](VectorXd MuBar_in){
    MatrixXd Hj = MatrixXd::Constant(3, 4, 0.0);

    float px = MuBar_in(0);
    float py = MuBar_in(1);
    float vx = MuBar_in(2);
    float vy = MuBar_in(3);

    float p_mag = px*px + py*py;
    float p_dot = sqrt(p_mag);

    if (p_mag < 0.001) {
      p_mag = 0.0001;  // arbitrary small number
      cout << "Warning: div0 in RadarEKF::Initialize encountered\n";
    }

    Hj << px/p_dot, py/p_dot, 0, 0,
         -py/p_mag, px/p_mag, 0, 0,
         py*(vx*py-vy*px)/(p_mag*p_dot), px*(vy*px-vx*py)/(p_mag*p_dot), px/p_dot, py/p_dot;

    return Hj;
  };

  MatrixXd Qt = MatrixXd::Constant(3, 3, 0.0);
  Qt << 0.09,0     ,0,
        0   ,0.0009,0,
        0   ,0     ,0.09;

  ExtendedKF::Initialize(g_in, Gt, Rt, h_in, H_in, Qt);
}
