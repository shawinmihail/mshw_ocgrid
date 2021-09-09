#include<iostream>
#include<mshw_ocgrid.h>
#include<cmath>
#include<vector>
#include <eigen3/Eigen/src/Core/ArrayBase.h>
#include"Eigen/Dense"
#include"trajectory_keeper.h"

int main(int argc, char *argv[])
{
    OcGrid ocGrid(10, 2, -50, -10, 10, 50);
    
    std::vector<Eigen::Vector2f> p1s;
    std::vector<Eigen::Vector2f> p2s;
    
    p1s.push_back(Eigen::Vector2f(0,0));
    //p1s.push_back(Eigen::Vector2f(0,0));
    //p1s.push_back(Eigen::Vector2f(0,0));
    //p2s.push_back(Eigen::Vector2f(1,0));
    //p2s.push_back(Eigen::Vector2f(1,1));
    p2s.push_back(Eigen::Vector2f(1,-1));
    
    float lim = 10; float free_addition = -1; float busy_addition = 5; 

    
    ocGrid.update(p1s, p2s, lim, free_addition, busy_addition, MapLayer::COMMON);
    
    std::vector<Eigen::Vector2f> vectors = ocGrid.get_obstcl_points(1, MapLayer::COMMON);
    for (auto v : vectors)
    {
        std::cout << std::endl << v.transpose() << std::endl;
    }

    ocGrid.print_map(COMMON);
    
    return 0;
}
