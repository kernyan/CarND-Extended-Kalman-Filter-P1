#ifndef UNSCENTED_KF_H_
#define UNSCENTED_KF_H_
#include "parametric_kf.h"
#include "measurement_package.h"
#include <functional>


class UnscentedKF : public ParametricKF{

public:

  UnscentedKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in);

  void Initialize(
      Func2 &g_in,
      MatrixXd &Rt_in,
      MatrixXd &nu_in,
      Func2 &h_in,
      MatrixXd &Qt_in);
  virtual void Step(MeasurementPackage &meas_in) override;

protected:

  void InitWeights();
  MatrixXd GenerateSigPts(const VectorXd &x_in,
      const MatrixXd &Sig_in) const;
  MatrixXd GenerateAugSigPts(const VectorXd &m_in,
      const MatrixXd &P_in,
      const MatrixXd &nu_in) const;
  MatrixXd Func_SigPtsSet(Func2 &f_in, double dt,
      const MatrixXd &SigSet_in) const;
  VectorXd PredictSigPts(double dt,
      const VectorXd &SigPts_in) const;
  VectorXd WeightedMean(const MatrixXd &SigPtsSet_in) const;
  MatrixXd CentralMoment(const MatrixXd &a_in,
      const VectorXd &mean_in) const;
  MatrixXd Covariance(const MatrixXd &a_in,
      const MatrixXd &b_in) const;
  void NormToPi(int Col_in, MatrixXd Mat_in) const;

  // member variables

  VectorXd &Mu_;    // shared with all filters
  MatrixXd &Sigma_; // shared with all filters
  long long &previous_time_; // shared with all filters

  int n_x_;       // state dimension
  int n_aug_;     // augmented state dimension
  double lambda_; // sigma point spreading parameter
  VectorXd weights_; // sigma point weights

  Func2 g_; // TransitionFunc
  MatrixXd Rt_; // TransitionSigma
  MatrixXd nu_; // process noise
  Func2 h_;      // MeasurementFunction
  MatrixXd Qt_; // MeasurementSigma
};


#endif /* UNSCENTED_KF_H_ */
