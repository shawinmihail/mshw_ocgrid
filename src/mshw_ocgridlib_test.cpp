#include<iostream>
#include<mshw_ocgrid.h>
#include<cmath>
#include<vector>
#include <eigen3/Eigen/src/Core/ArrayBase.h>
#include"Eigen/Dense"
#include"trajectory_keeper.h"

int main(int argc, char *argv[])
{
    OcGrid ocGrid(15, 2, -50, -10, 10, 50);
    ocGrid.fill_for_test(COMMON);
    ocGrid.print_map(COMMON);
    
    auto local_map = ocGrid.get_local_map(Eigen::Vector2f(4,0), 1.1, COMMON);
    std::cout << std::endl << local_map.map.transpose().cast<int>() << std::endl;
    std::cout << std::endl << local_map.size << std::endl;
    std::cout << std::endl << local_map.resolution << std::endl;
    std::cout << std::endl << local_map.dimention << std::endl;
    

    return 0;
}
