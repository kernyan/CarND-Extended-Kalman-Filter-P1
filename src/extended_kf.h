#ifndef EXTENDED_KF_H_
#define EXTENDED_KF_H_
#include "parametric_kf.h"
#include "measurement_package.h"
#include <functional>


class ExtendedKF : public ParametricKF {

public:

  ExtendedKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in);

  void Initialize(Func &g_in,
      MatrixXd &G_in,
      MatrixXd &Rt_in,
      Func &h_in,
      function<MatrixXd (VectorXd)> &H_in,
      MatrixXd &Qt_in);
  virtual void Step(MeasurementPackage &meas_in) override;
  VectorXd GetMuBar() const {return MuBar_;};
  MatrixXd GetSigmaBar() const {return SigmaBar_;};


protected:

  virtual void CalculateMuBar();
  virtual void CalculateSigmaBar();
  virtual MatrixXd CalculateMeasurementVar();
  virtual VectorXd CalculatePredictedMeasurement();
  MatrixXd CalculatePredJacobian() const;
  MatrixXd CalculateMeasJacobian() const;

  // member variables

  VectorXd &Mu_;    // shared with all filters
  MatrixXd &Sigma_; // shared with all filters
  long long &previous_time_; // shared with all filters

  VectorXd MuBar_;
  MatrixXd SigmaBar_;
  Func g_;      // TransitionFunction
  MatrixXd G_;  // Jacobian of TransitionFunction
  MatrixXd Rt_; // TransitionSigma
  Func h_;      // MeasurementFunction
  function<MatrixXd (VectorXd)> H_;  // Jacobian of MeasurementFunction
  MatrixXd Qt_; // MeasurementSigma
};


#endif /* EXTENDED_KF_H_ */
