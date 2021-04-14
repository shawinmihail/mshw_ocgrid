#include<iostream>
#include<mshw_ocgrid.h>
#include<cmath>
#include<vector>
#include"Eigen/Dense"

int main(int argc, char *argv[])
{
    OcGrid og(50, 0.5);
    

    std::vector<Eigen::Vector2f> p1s;
    p1s.push_back(Eigen::Vector2f(0.7, 0.0));
    std::vector<Eigen::Vector2f> p2s;
    p2s.push_back(Eigen::Vector2f(6.04, 0.93));
    
    og.update(p1s, p2s);
    std::vector<Eigen::Vector2f> vectors = og.get_obstcl_poits();
    for (auto v : vectors)
    {
        std::cout << std::endl << v.transpose() << std::endl;
    }

    //std::cout << std::endl << og.get_map().cast<int>() << std::endl;
    
    return 0;
}
