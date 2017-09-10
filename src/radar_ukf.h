#ifndef RADAR_UKF_H_
#define RADAR_UKF_H_
#include "unscented_kf.h"
#include "measurement_package.h"
#include <functional>


class RadarUKF : public UnscentedKF {

  public:

    RadarUKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in);
    virtual void Step(MeasurementPackage &meas_in) override final;
    void Initialize();

  private:

    // member variables

    const SensorType Sensor_ = RADAR;
};


#endif /* RADAR_UKF_H_ */
