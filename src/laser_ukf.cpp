#include "laser_ukf.h"
#include <iostream>


void LaserUKF::Step(MeasurementPackage &meas_in){

  if (Sensor_ != meas_in.sensor_type_) return;

  //UnscentedKF::Step(meas_in);

  double dt = meas_in.timestamp_ - previous_time_;
  previous_time_ = meas_in.timestamp_;

  dt *= 0.000001; // scale to seconds

  MatrixXd AugSigPts = GenerateAugSigPts(Mu_, Sigma_, nu_);
  MatrixXd XSig_Pred = Func_SigPtsSet(g_, dt, AugSigPts);
  VectorXd Mu_Bar = WeightedMean(XSig_Pred);
  MatrixXd FirstMoment = CentralMoment(XSig_Pred, Mu_Bar);
  NormToPi(3, FirstMoment); // yaw to within +/- pi
  MatrixXd Sigma_Bar=Covariance(FirstMoment,FirstMoment);

  // Linear measurement update since laser meas is linear
  MatrixXd H = MatrixXd::Identity(2,5);
  VectorXd Zhat = H * Mu_Bar;
  VectorXd y = meas_in.raw_measurements_ - Zhat;
  MatrixXd S = H * Sigma_Bar * H.transpose() + Qt_;
  MatrixXd K = Sigma_Bar * H.transpose() * S.inverse();

  Mu_ = Mu_Bar + K * y;
  Sigma_ = Sigma_Bar - K * S * K.transpose();
}


void LaserUKF::Initialize(){

  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_aug_;

  InitWeights();

  auto n_x = n_x_;
  Func2 g_in = [n_x](double dt, const VectorXd &x_in){
    double px       = x_in(0);
    double py       = x_in(1);
    double v        = x_in(2);
    double yaw      = x_in(3);
    double yawd     = x_in(4);
    double nu_a     = x_in(5);
    double nu_yawdd = x_in(6);

    double px_p, py_p, v_p, yaw_p, yawd_p = 0;

    if (fabs(yawd) > 0.001){
      px_p = px+(v/yawd)*( sin(yaw+yawd*dt) - sin(yaw));
      py_p = py+(v/yawd)*(-cos(yaw+yawd*dt) + cos(yaw));
    } else {
      px_p = v * dt * cos(yaw);
      py_p = v * dt * sin(yaw);
    }

    px_p  += 0.5 * dt * dt * nu_a * cos(yaw);
    py_p  += 0.5 * dt * dt * nu_a * sin(yaw);
    v_p    = v + dt * nu_a;
    yaw_p  = yaw + yawd * dt + 0.5 * dt * dt * nu_yawdd;
    yawd_p = yawd + dt * nu_yawdd;

    VectorXd Out = VectorXd::Constant(n_x, 0.0);
    Out << px_p, py_p, v_p, yaw_p, yawd_p;

    return Out;
  };

  MatrixXd Rt_in = MatrixXd::Identity(n_x_, n_x_);
  Rt_in(0,0) = 0.15;
  Rt_in(1,1) = 0.15;

  double noise_a = 4;
  double noise_yawdd = 0.09;
  int Dim = n_aug_ - n_x_;
  MatrixXd nu_in = MatrixXd::Constant(Dim, Dim, 0.0);
  nu_in(0,0) = noise_a; // longitud acceleration process var
  nu_in(1,1) = noise_yawdd; // yaw acceleration process var

  Func2 h_in; // not used

  MatrixXd Qt_in(2,2);
  Qt_in << 0.0225, 0     ,
           0     , 0.0225;


  UnscentedKF::Initialize(g_in, Rt_in, nu_in, h_in, Qt_in);
}


LaserUKF::LaserUKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in) :
  UnscentedKF(Mu_in, Sigma_in, t_in)
{
}
