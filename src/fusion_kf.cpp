#include "fusion_kf.h"
#include "Eigen/Dense"
#include <iostream>
#include <vector>


FusionKF::FusionKF(SensorModel Model_in) :
  IsFirstTime (true),
  previous_time_(-1),
  KinematicModel (Model_in)
{
  switch (Model_in){

  case CONSTANT_VELOCITY:

    Mu_ = VectorXd(4);
    Mu_ << 0,0,0,0;

    Sigma_ = MatrixXd(4,4);
    Sigma_ << 1,0,0,0,
              0,1,0,0,
              0,0,1000,0,
              0,0,0,1000;

    break;

  case CONSTANT_TURNRATE_VELOCITY:

    Mu_ = VectorXd(5);
    Mu_ << 1,1,1,1,0.1;

    Sigma_ = MatrixXd(5,5);
    Sigma_ << 0.15,0   ,0,0,0,
              0   ,0.15,0,0,0,
              0   ,0   ,1,0,0,
              0   ,0   ,0,1,0,
              0   ,0   ,0,0,1;

    break;

  default:

    assert(0);
  }
}


FusionKF::~FusionKF(){
  for (auto &each : filters_){
    delete each;
    each = nullptr;
  }
}

void FusionKF::AddLaserLKF(){
  
  // Laser Filter
  assert(KinematicModel == CONSTANT_VELOCITY); // LaseLKF only supported with Constant Velocity model

  LaserLKF* LaserFilter = new LaserLKF(Mu_, Sigma_, previous_time_);
  LaserFilter->Initialize();
  filters_.push_back(LaserFilter);
}

void FusionKF::AddRadarEKF(){
  
  // Radar Filter
  
  assert(KinematicModel == CONSTANT_VELOCITY); // RadarEKF only supported with Constant Velocity model

  RadarEKF *RadarFilter = new RadarEKF(Mu_, Sigma_, previous_time_);
  RadarFilter->Initialize();
  filters_.push_back(RadarFilter);
}


void FusionKF::AddLaserUKF(){
  
  // Radar Filter
  
  assert(KinematicModel == CONSTANT_TURNRATE_VELOCITY); // LaserUKF only supported with CTRV model

  LaserUKF *LaserFilter = new LaserUKF(Mu_, Sigma_, previous_time_);
  LaserFilter->Initialize();
  filters_.push_back(LaserFilter);
}


void FusionKF::AddRadarUKF(){
  
  // Radar Filter
  
  assert(KinematicModel == CONSTANT_TURNRATE_VELOCITY); // RadarUKF only supported with CTRV model

  RadarUKF *RadarFilter = new RadarUKF(Mu_, Sigma_, previous_time_);
  RadarFilter->Initialize();
  filters_.push_back(RadarFilter);
}


void FusionKF::ProcessMeasurement(MeasurementPackage &meas_in){
  if (IsFirstTime){
    Mu_(0) = meas_in.raw_measurements_(0);
    Mu_(1) = meas_in.raw_measurements_(1);
    assert(previous_time_ == -1);
    previous_time_ = meas_in.timestamp_;
    
    IsFirstTime = false;
  }

  for (auto &filter : filters_){
    filter->Step(meas_in);
  }
}


VectorXd CalculateRMSE(const vector<VectorXd> &estimations,
                       const vector<VectorXd> &ground_truth){
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if (estimations.size() != ground_truth.size()
     || estimations.size() == 0){
    cout << "Invalud estimation or ground truth data" << endl;
    return rmse;
  }

  for (int i = 0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

