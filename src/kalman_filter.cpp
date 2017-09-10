#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

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

  ExtendedKF::Step(meas_in);
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
    if (phi > -M_PI){
      while (phi > M_PI) phi -= 2*M_PI;
    } else {
      while (phi < -M_PI) phi += 2*M_PI;
    }

    float p_dot = (px*vx+py*vy)/p;
    zhat << p, phi, p_dot;
    return zhat;
  };

  function<MatrixXd (VectorXd)> H_in = [](VectorXd MuBar_in){
    MatrixXd Hj(3,4);

    float px = MuBar_in(0);
    float py = MuBar_in(1);
    float vx = MuBar_in(2);
    float vy = MuBar_in(3);

    float p_mag = px*px + py*py;
    float p_dot = sqrt(p_mag);

    if (p_mag < 0.001) return Hj;

    Hj << px/p_dot, py/p_dot, 0, 0,
         -py/p_mag, px/p_mag, 0, 0,
         py*(vx*py-vy*px)/(p_mag*p_dot), px*(vy*px-vx*py)/(p_mag*p_dot), px/p_dot, py/p_dot;
    
    return Hj;
  };

  MatrixXd Qt(3,3);
  Qt << 0.09,0     ,0,
        0   ,0.0009,0,
        0   ,0     ,0.09;

  ExtendedKF::Initialize(g_in, Gt, Rt, h_in, H_in, Qt);
}


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
