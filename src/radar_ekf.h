#ifndef RADAR_EKF_H_
#define RADAR_EKF_H_
#include "Eigen/Dense"
#include "extended_kf.h"
#include "measurement_package.h"
#include <functional>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class RadarEKF : public ExtendedKF {

  public:

    RadarEKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in);
    virtual void Step(MeasurementPackage &meas_in) override final;
    void Initialize();

  private:

    const SensorType Sensor_ = RADAR;
    float noise_ax_;
    float noise_ay_;
};

#endif /* RADAR_EKF_H_ */
