#ifndef LASER_UKF_H_
#define LASER_UKF_H_
#include "unscented_kf.h"
#include "measurement_package.h"
#include <functional>


class LaserUKF : public UnscentedKF {

  public:

    LaserUKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in);
    virtual void Step(MeasurementPackage &meas_in) override final;
    void Initialize();

  private:

    // member variables

    const SensorType Sensor_ = LASER;
};

#endif /* LASER_UKF_H_ */
