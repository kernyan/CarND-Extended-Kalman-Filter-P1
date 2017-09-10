#ifndef LINEAR_KF_H_
#define LINEAR_KF_H_
#include "parametric_kf.h"
#include "measurement_package.h"
#include <functional>

class LinearKF : public ParametricKF {

public:

  LinearKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in);
  void Initialize(MatrixXd &TransFunc_in,
      MatrixXd &TransSigma_in,
      MatrixXd &MeasFunc_in,
      MatrixXd &MeasSigma_in);
  virtual void Step(MeasurementPackage &meas_in) override;
  VectorXd GetMuBar() const {return MuBar_;};
  MatrixXd GetSigmaBar() const {return SigmaBar_;};

protected:

  virtual void CalculateMuBar();
  virtual void CalculateSigmaBar();
  virtual MatrixXd CalculateMeasurementVar(); // S
  virtual VectorXd CalculatePredictedMeasurement(); // zhat

  // member variables

  VectorXd &Mu_;    // shared with all filters
  MatrixXd &Sigma_; // shared with all filters
  long long &previous_time_; // shared with all filters

  VectorXd MuBar_;
  MatrixXd SigmaBar_;
  MatrixXd At_; // TransitionFunction
  MatrixXd Rt_; // TransitionSigma
  MatrixXd Ct_; // MeasurementFunction
  MatrixXd Qt_; // MeasurementSigma
};


#endif /* LINEAR_KF_H_ */
