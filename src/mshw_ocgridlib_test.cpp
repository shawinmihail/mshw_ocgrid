#include<iostream>
#include<mshw_ocgrid.h>
#include<cmath>
#include<vector>
#include"Eigen/Dense"

int main(int argc, char *argv[])
{
    OcGrid og(10, 0.1);
    

    std::vector<Eigen::Vector2f> p1s;
    p1s.push_back(Eigen::Vector2f(0,0));
    p1s.push_back(Eigen::Vector2f(0,0));
    p1s.push_back(Eigen::Vector2f(0,0));
    std::vector<Eigen::Vector2f> p2s;
    p2s.push_back(Eigen::Vector2f(1,2));
    p2s.push_back(Eigen::Vector2f(2,1));
    p2s.push_back(Eigen::Vector2f(-1.75,2.25));
    
    og.update(p1s, p2s);

    std::cout << std::endl << og.get_map().cast<int>() << std::endl;
    
    return 0;
}
