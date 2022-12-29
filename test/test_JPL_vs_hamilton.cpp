
#include "extern/JPL_imu_error.h"
#include "extern/vinsmono_imu_error.h"
#include "PoseSpline/maplab/maplab_imu_factor.h"
#include "internal/pose_local_parameterization.h"
#include "PoseSpline/PoseLocalParameter.hpp"
#include "PoseSpline/NumbDifferentiator.hpp"
#include "PoseSpline/maplab/pose_param_jpl.h"


#include "PoseSpline/Pose.hpp"
#include "PoseSpline/PoseSpline.hpp"
#include "PoseSpline/Time.hpp"
#include "csv.h"

#include <ceres/gradient_checker.h>


int main() {

    Eigen::Vector3d un_gyr(0.1, 0.32, 0.41);
    double dt = 1;
    Eigen::Quaterniond hamilton_q_WI0(0.3, 0.9, -0.4, 0.5);
    hamilton_q_WI0.normalize();
    std::cout << "hamilton_q_WI0: " << hamilton_q_WI0.coeffs().transpose() << std::endl;
    Eigen::Quaterniond hamilton_q_WI1 =
                    hamilton_q_WI0 * Eigen::Quaterniond(1, un_gyr(0) * dt / 2, un_gyr(1) * dt / 2, un_gyr(2) * dt / 2);
    hamilton_q_WI1.normalize();
    std::cout << "hamilton_q_WI1: " << hamilton_q_WI1.coeffs().transpose() << std::endl;

    QuaternionTemplate<double> JPL_q_WI0;
    common::fromRotationMatrixJPL(hamilton_q_WI0.toRotationMatrix(), &JPL_q_WI0);
    QuaternionTemplate<double> JPL_q_I0W = common::quaternionInverseJPL(JPL_q_WI0);
    QuaternionTemplate<double> JPL_q_I1W;
    Eigen::Vector3d delta(un_gyr * dt);
    // delta = -delta;
    common::positiveQuaternionProductJPL(deltaQuat<double>(delta),JPL_q_I0W, JPL_q_I1W);
    JPL_q_I1W = JPL_q_I1W/JPL_q_I1W.norm();

    std::cout << "JPL_q_I0W: " << JPL_q_I0W.transpose() << std::endl;
    std::cout << "JPL_q_I1W: " << JPL_q_I1W.transpose() << std::endl;



    return 0;

}