#ifndef FusionKF_H_
#define FusionKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include "laser_lkf.h"
#include "radar_ekf.h"
#include "radar_ukf.h"
#include "laser_ukf.h"

using namespace std;

enum SensorModel{
  CONSTANT_VELOCITY,
  CONSTANT_TURNRATE_VELOCITY
};


class FusionKF {
public:
  
  list<ParametricKF*> filters_;
  void AddLaserLKF();
  void AddRadarEKF();
  void AddLaserUKF();
  void AddRadarUKF();
  void ProcessMeasurement(MeasurementPackage &meas_in);
  VectorXd GetMu() const {return Mu_;};
  MatrixXd GetSigma() const {return Sigma_;};

  FusionKF() = delete;
  FusionKF(SensorModel Model_in);
  ~FusionKF();
  const SensorModel KinematicModel;

private:

  VectorXd Mu_;
  MatrixXd Sigma_;
  long long previous_time_;
  bool IsFirstTime;
};

VectorXd CalculateRMSE(const vector<VectorXd> &estimations,
                       const vector<VectorXd> &ground_truth);
#endif /* FUSION_KF_H */
