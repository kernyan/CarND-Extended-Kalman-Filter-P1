#include "unscented_kf.h"
#include <iostream>


UnscentedKF::UnscentedKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in) :
    Mu_ (Mu_in),
    Sigma_ (Sigma_in),
    previous_time_ (t_in)
{
}


void UnscentedKF::Initialize(
    Func2 &g_in,
    MatrixXd &Rt_in,
    MatrixXd &nu_in,
    Func2 &h_in,
    MatrixXd &Qt_in){
  g_ = g_in;
  Rt_ = Rt_in;
  nu_ = nu_in;
  h_ = h_in;
  Qt_ = Qt_in;
}


void UnscentedKF::InitWeights(){

  double w_i = 0.5 / (n_aug_ + lambda_);
  weights_ = VectorXd(1 + 2*n_aug_);
  weights_.fill(w_i);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
}


MatrixXd UnscentedKF::GenerateAugSigPts(const VectorXd &m_in,
    const MatrixXd &P_in,
    const MatrixXd &nu_in) const {
  VectorXd x_aug = VectorXd::Constant(n_aug_, 0.0);
  x_aug.head(n_x_) = m_in; // assumes 0 mean of noise process

  int kAugSigPts = 1 + n_aug_ * 2;
  MatrixXd P_aug = MatrixXd::Constant(n_aug_, n_aug_, 0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_in;
  P_aug.bottomRightCorner(n_aug_-n_x_, n_aug_-n_x_) = nu_in;
  return GenerateSigPts(x_aug, P_aug);
}


MatrixXd UnscentedKF::GenerateSigPts(const VectorXd &x_in,
    const MatrixXd &Sig_in) const {
  int Dim = x_in.size();
  int kSigPts = 1 + Dim * 2; // 1 for mean, 2 for each +/-

  MatrixXd Out = MatrixXd::Constant(Dim, kSigPts, 0.0);
  MatrixXd A = Sig_in.llt().matrixL();

  double Dispersion = sqrt(lambda_ + Dim);
  MatrixXd gamma = A * Dispersion;

  Out.col(0) = x_in;
  for (int i = 0; i < Dim; ++i){
    Out.col(i+1)     = x_in + gamma.col(i);
    Out.col(i+1+Dim) = x_in - gamma.col(i);
  }

  return Out;
}


MatrixXd UnscentedKF::Func_SigPtsSet(Func2 &f_in, double dt,
    const MatrixXd &SigSet_in) const {

  int kSigPtsSet = SigSet_in.cols();

  VectorXd Temp = f_in(dt, SigSet_in.col(0));
  int Dim = Temp.rows();

  MatrixXd Out = MatrixXd::Constant(Dim, kSigPtsSet, 0.0);

  for (int i = 0; i < kSigPtsSet; ++i){
    Out.col(i) = f_in(dt, SigSet_in.col(i));
  }

  return Out;
}


VectorXd UnscentedKF::PredictSigPts(double dt,
    const VectorXd &SigPts_in) const {

  return g_(dt, SigPts_in);
}


VectorXd UnscentedKF::WeightedMean(
    const MatrixXd &SigPtsSet_in) const {

  int kState = SigPtsSet_in.rows();
  int kSigPtsSet = SigPtsSet_in.cols();

  VectorXd Out = VectorXd::Constant(kState, 0.0);
  for (int i = 0; i < kSigPtsSet; ++i){
    Out += weights_(i) * SigPtsSet_in.col(i);
  }

  return Out;
}


MatrixXd UnscentedKF::CentralMoment(const MatrixXd &a_in,
    const VectorXd &mean_in) const {

  int kState = mean_in.rows();
  int kSigPtsSet = a_in.cols();

  MatrixXd Out = MatrixXd::Constant(kState, kSigPtsSet, 0.0);
  for (int i = 0; i < kSigPtsSet; ++i){
    Out.col(i) = a_in.col(i) - mean_in;
  }

  return Out;
}


MatrixXd UnscentedKF::Covariance(const MatrixXd &a_in,
    const MatrixXd &b_in) const {

  int kDimIn = a_in.rows();
  int kDimOut = b_in.rows();
  int kSigPtsSet = a_in.cols();
  assert(kSigPtsSet == b_in.cols());

  MatrixXd Out = MatrixXd::Constant(kDimIn, kDimOut, 0.0);
  for (int i = 0; i < kSigPtsSet; ++i){
    Out += weights_(i) * a_in.col(i)*b_in.col(i).transpose();
  }

  return Out;
}


void UnscentedKF::NormToPi(int Row_in, MatrixXd Mat_in) const {

  // normalize value to within pi range

 int Dim = Mat_in.cols();

 for (int i = 0; i < Dim; ++i){
   while (Mat_in.col(i)(Row_in) > M_PI)
     Mat_in.col(i)(Row_in) -= 2.*M_PI;
   while (Mat_in.col(i)(Row_in) < -M_PI)
     Mat_in.col(i)(Row_in) += 2.*M_PI;
 }
}


void UnscentedKF::Step(MeasurementPackage &meas_in){

  double dt = meas_in.timestamp_ - previous_time_;
  previous_time_ = meas_in.timestamp_;

  dt *= 0.000001; // scale to seconds

  cout << "dt " << dt << endl;

  MatrixXd AugSigPts = GenerateAugSigPts(Mu_, Sigma_, nu_);
  MatrixXd XSig_Pred = Func_SigPtsSet(g_, dt, AugSigPts);
  VectorXd Mu_Bar = WeightedMean(XSig_Pred);
  MatrixXd FirstMoment = CentralMoment(XSig_Pred, Mu_Bar);
  MatrixXd Sigma_Bar=Covariance(FirstMoment,FirstMoment);

  // common to replace update step with other implementation
  // e.g. 1. with reusing transition prediction sigma pts, or
  //      2. with Linear measurement update if linear
  MatrixXd nu_0 = MatrixXd::Zero(2,2);
  MatrixXd SigPts = GenerateAugSigPts(Mu_Bar, Sigma_Bar,nu_0);
  MatrixXd ZSig_Pred = Func_SigPtsSet(h_, dt, SigPts);
  VectorXd Zhat = WeightedMean(ZSig_Pred);
  MatrixXd FirstMoment2 = CentralMoment(ZSig_Pred, Zhat);
  MatrixXd S = Covariance(FirstMoment2, FirstMoment2)+Qt_;
  MatrixXd CrossCov = Covariance(FirstMoment, FirstMoment2);
  MatrixXd K = CrossCov * S.inverse();
  VectorXd y = meas_in.raw_measurements_ - Zhat;

  Mu_ = Mu_Bar + K * y;
  Sigma_ = Sigma_Bar - K * S * K.transpose();
}


