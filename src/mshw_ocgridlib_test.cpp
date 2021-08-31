#include<iostream>
#include<mshw_ocgrid.h>
#include<cmath>
#include<vector>
#include <eigen3/Eigen/src/Core/ArrayBase.h>
#include"Eigen/Dense"
#include"trajectory_keeper.h"

int main(int argc, char *argv[])
{
    TrajectoryKeeper tk;
    tk.reinit(3);
    
    float t0 = 0;
    Eigen::Vector3f p0(0.0, 0.0, 0.0);
    Eigen::Vector4f q0(1.0, 0.0, 0.0, 0.0);
    for (int i = 0; i < 7; i++)
    {
        std::cout << i << std::endl;
        t0 = t0 + 0.1;
        p0 = p0 + Eigen::Vector3f(1,2,3);
        q0 = q0 + Eigen::Vector4f(0,0,0,0.1);
        q0 = q0 / q0.norm();
        bool full = tk.addPoint(t0, p0, q0);
        if (full)
        {
            std::cout << "full" << std::endl;
            break;
        }
    }
    
    TrajectoryKeeperPoint p_out;
    bool res = tk.getApproximation(0.3, p_out);
    if(res)
        std::cout << "OK" << std::endl;
    else
        std::cout << "NOT OK" << std::endl;
    p_out.print();
}
