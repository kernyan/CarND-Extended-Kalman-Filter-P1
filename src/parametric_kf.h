#ifndef PARAMETRIC_KF_H_
#define PARAMETRIC_KF_H_  
#include "Eigen/Dense"
#include "measurement_package.h"
#include <functional>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using Func = function<VectorXd (VectorXd)>;
using Func2 = function<VectorXd (double, const VectorXd &)>;

class ParametricKF {

  // abstract interface for all other filter implementation

public:

	virtual ~ParametricKF(){};
	virtual void Step(MeasurementPackage &meas_in) = 0;
};



#endif /* PARAMETRIC_KF_H_ */

