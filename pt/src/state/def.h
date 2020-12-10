#pragma once
#include <Eigen/Core>

using ScalarList = Eigen::ArrayXd;
using VectorList = Eigen::ArrayX3d;
using Vector = Eigen::RowVector3d;
#define iter(vl) vl.rowwise()
