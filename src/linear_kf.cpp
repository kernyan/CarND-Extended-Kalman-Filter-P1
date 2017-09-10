#include "linear_kf.h"
#include <iostream>


LinearKF::LinearKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in) :
  Mu_ (Mu_in),
  Sigma_ (Sigma_in),
  previous_time_(t_in)
{
}


void LinearKF::Initialize(MatrixXd &TransFunc_in,
  MatrixXd &TransSigma_in,
  MatrixXd &MeasFunc_in,
  MatrixXd &MeasSigma_in){

  At_ = TransFunc_in;
  Rt_ = TransSigma_in;
  Ct_ = MeasFunc_in;
  Qt_ = MeasSigma_in;
}


void LinearKF::CalculateMuBar(){
  MuBar_ = At_ * Mu_;
}


void LinearKF::CalculateSigmaBar(){
  SigmaBar_ = At_ * Sigma_ * At_.transpose() + Rt_;
}


MatrixXd LinearKF::CalculateMeasurementVar(){
  return Ct_ * SigmaBar_ * Ct_.transpose() + Qt_;
}


VectorXd LinearKF::CalculatePredictedMeasurement(){
  return Ct_ * MuBar_;
}


void LinearKF::Step(MeasurementPackage &meas_in){
  CalculateMuBar();
  CalculateSigmaBar();
  MatrixXd S = CalculateMeasurementVar();
  MatrixXd K = SigmaBar_ * Ct_.transpose() *  S.inverse();
  VectorXd zhat = CalculatePredictedMeasurement();
  Mu_ = MuBar_ + K * (meas_in.raw_measurements_ - zhat);
  Sigma_ = SigmaBar_ - K * S * K.transpose();
}

