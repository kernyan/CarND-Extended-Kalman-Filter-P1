#include "tools.h"


void NormToPi(int Row_in, MatrixXd Mat_in) {

  // normalize value to within pi range

 int Dim = Mat_in.cols();

 for (int i = 0; i < Dim; ++i){
   while (Mat_in.col(i)(Row_in) > M_PI)
     Mat_in.col(i)(Row_in) -= 2.*M_PI;
   while (Mat_in.col(i)(Row_in) < -M_PI)
     Mat_in.col(i)(Row_in) += 2.*M_PI;
 }
}
