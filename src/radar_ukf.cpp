#include "radar_ukf.h"
#include <iostream>


void RadarUKF::Step(MeasurementPackage &meas_in){

  if (Sensor_ != meas_in.sensor_type_) return;

  double dt = meas_in.timestamp_ - previous_time_;
  previous_time_ = meas_in.timestamp_;

  dt *= 0.000001; // scale to seconds


  MatrixXd AugSigPts = GenerateAugSigPts(Mu_, Sigma_, nu_);
  MatrixXd XSig_Pred = Func_SigPtsSet(g_, dt, AugSigPts);
  VectorXd Mu_Bar = WeightedMean(XSig_Pred);
  MatrixXd FirstMoment = CentralMoment(XSig_Pred, Mu_Bar);
  NormToPi(3, FirstMoment); // yaw to within +/- pi
  MatrixXd Sigma_Bar=Covariance(FirstMoment,FirstMoment);

  // Reusing transition prediction sigma pts
  MatrixXd ZSig_Pred = Func_SigPtsSet(h_, dt, XSig_Pred);
  VectorXd Zhat = WeightedMean(ZSig_Pred);
  MatrixXd FirstMoment2 = CentralMoment(ZSig_Pred, Zhat);
  NormToPi(1, FirstMoment2); // phi to within +/- pi
  MatrixXd S = Covariance(FirstMoment2, FirstMoment2)+Qt_;
  MatrixXd CrossCov = Covariance(FirstMoment, FirstMoment2);
  MatrixXd K = CrossCov * S.inverse();
  VectorXd y = meas_in.raw_measurements_ - Zhat;

  while (y(1) >  M_PI) y(1) -= 2.*M_PI;
  while (y(1) < -M_PI) y(1) += 2.*M_PI;

  Mu_ = Mu_Bar + K * y;
  Sigma_ = Sigma_Bar - K * S * K.transpose();
}


void RadarUKF::Initialize(){

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

  Func2 h_in = [](double dt, const VectorXd x_in){

    short kMeasurement = 3; // for rho, phi, and rho dot
    double px   = x_in(0);
    double py   = x_in(1);
    double v    = x_in(2);
    double yaw  = x_in(3);
    double vx   = v * cos(yaw);
    double vy   = v * sin(yaw);

    VectorXd Out = VectorXd::Constant(kMeasurement, 0.0);

    Out(0) = sqrt(px*px + py*py);
    Out(1) = atan2(py, px);
    Out(2) = (px * vx + py * vy)/Out(0);

    return Out;
  };

  MatrixXd Qt_in(3,3);
  Qt_in << 0.09,0     ,0,
           0   ,0.0009,0,
           0   ,0     ,0.09;

  UnscentedKF::Initialize(g_in, Rt_in, nu_in, h_in, Qt_in);
}


RadarUKF::RadarUKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in) :
  UnscentedKF(Mu_in, Sigma_in, t_in)
{
}

