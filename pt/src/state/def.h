#pragma once
#include <Eigen/Core>

using ScalarArray = Eigen::ArrayXd;
using Vec3DArray = Eigen::ArrayX3d;
using Vec3D = Eigen::RowVector3d;
#define iter(vl) vl.rowwise()
