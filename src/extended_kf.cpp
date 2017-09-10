#include "extended_kf.h"
#include <iostream>


ExtendedKF::ExtendedKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in) :
    Mu_ (Mu_in),
    Sigma_ (Sigma_in),
    previous_time_ (t_in)
{
}


void ExtendedKF::Initialize(Func &g_in,
    MatrixXd &G_in,
    MatrixXd &Rt_in,
    Func &h_in,
    function<MatrixXd (VectorXd)> &H_in,
    MatrixXd &Qt_in){
  g_ = g_in;
  G_ = G_in;
  Rt_ = Rt_in;
  h_ = h_in;
  H_ = H_in;
  Qt_ = Qt_in;
}


void ExtendedKF::CalculateMuBar(){
  MuBar_ = g_(Mu_);
}


void ExtendedKF::CalculateSigmaBar(){
  SigmaBar_ = G_ * Sigma_ * G_.transpose() + Rt_;
}


MatrixXd ExtendedKF::CalculateMeasurementVar(){
  MatrixXd Ht = H_(MuBar_);
  return Ht * SigmaBar_ * Ht.transpose() + Qt_;
}

VectorXd ExtendedKF::CalculatePredictedMeasurement(){
  return h_(MuBar_);
}


void ExtendedKF::Step(MeasurementPackage &meas_in){
  CalculateMuBar();
  CalculateSigmaBar();
  MatrixXd S = CalculateMeasurementVar();
  MatrixXd Ht = H_(MuBar_);
  MatrixXd K = SigmaBar_ * Ht.transpose() * S.inverse();
  VectorXd zhat = CalculatePredictedMeasurement();
  VectorXd y = meas_in.raw_measurements_ - zhat;
  if (meas_in.raw_measurements_(1) * zhat(1) <= 0){
    y << 0,0,0; // do not update when theta flips sign
  }
  Mu_ = MuBar_ + K * y;
  Sigma_ = SigmaBar_ - K * S * K.transpose();
}


