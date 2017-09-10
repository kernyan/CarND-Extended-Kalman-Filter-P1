#ifndef LASER_LKF_H_
#define LASER_LKF_H_
#include "linear_kf.h"
#include "measurement_package.h"
#include <functional>


class LaserLKF : public LinearKF {

  public:

    LaserLKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in);
    virtual void Step(MeasurementPackage &meas_in) override final;
    void Initialize();

  private:

    const SensorType Sensor_ = LASER;
    float noise_ax_;
    float noise_ay_;
};


#endif /* LASER_LKF_H_ */
