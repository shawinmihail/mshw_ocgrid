#include<iostream>
#include<mshw_ocgrid.h>
#include<cmath>
#include<vector>
#include <eigen3/Eigen/src/Core/ArrayBase.h>
#include"Eigen/Dense"

int main(int argc, char *argv[])
{
    OcGrid og(3.5, 0.5);
    

    std::vector<Eigen::Vector2f> p1s;
    p1s.push_back(Eigen::Vector2f(0.0, 0.0));
    p1s.push_back(Eigen::Vector2f(0.0, 0.0));
    p1s.push_back(Eigen::Vector2f(0.0, 0.0));

    std::vector<Eigen::Vector2f> p2s;
    p2s.push_back(Eigen::Vector2f(1.1, 1.0));
    p2s.push_back(Eigen::Vector2f(1.1, 0.0));
    p2s.push_back(Eigen::Vector2f(1.25, 1.25));

    og.update(p1s, p2s, 10);
    std::vector<Eigen::Vector2f> vectors = og.get_obstcl_points();
    for (auto v : vectors)
    {
        std::cout << std::endl << v.transpose() << std::endl;
    }

    std::cout << std::endl << og.get_map().cast<int>() << std::endl;
    
    return 0;
}
