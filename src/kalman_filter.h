#ifndef PARAMETRIC_KF_H_
#define PARAMETRIC_KF_H_  
#include "Eigen/Dense"
#include "measurement_package.h"
#include <functional>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class ParametricKF {

  // abstract interface for all other filter implementation

public:

	virtual ~ParametricKF(){};
	virtual void Step(MeasurementPackage &meas_in) = 0;
};


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


using Func = function<VectorXd (VectorXd)>;

class ExtendedKF : public ParametricKF{
  
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


using Func2 = function<VectorXd (double, const VectorXd &)>;

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


class RadarUKF : public UnscentedKF {

  public:
    
    RadarUKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in);
    virtual void Step(MeasurementPackage &meas_in) override final;
    void Initialize();

  private:

    // member variables
  
    const SensorType Sensor_ = RADAR;
};


class LaserUKF : public UnscentedKF {

  public:
    
    LaserUKF(VectorXd &Mu_in, MatrixXd &Sigma_in, long long &t_in);
    virtual void Step(MeasurementPackage &meas_in) override final;
    void Initialize();

  private:

    // member variables
  
    const SensorType Sensor_ = LASER;
};
#endif /* PARAMETRIC_KF_H_ */

