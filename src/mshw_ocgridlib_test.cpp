#include<iostream>
#include<mshw_ocgrid.h>
#include<cmath>
#include<vector>
#include<set>
#include <eigen3/Eigen/src/Core/ArrayBase.h>
#include"Eigen/Dense"
#include"trajectory_keeper.h"

int main(int argc, char *argv[])
{
    /*
    OcGrid ocGrid(15, 2, -50, -10, 10, 50);
    ocGrid.fill_for_test(COMMON);
    ocGrid.print_map(COMMON);
    
    auto local_map = ocGrid.get_local_map(Eigen::Vector2f(4,0), 1.1, COMMON);
    std::cout << std::endl << local_map.map.transpose().cast<int>() << std::endl;
    std::cout << std::endl << local_map.size << std::endl;
    std::cout << std::endl << local_map.resolution << std::endl;
    std::cout << std::endl << local_map.dimention << std::endl;
    */
    
    std::vector<GridIndexWithVal> gs;
    gs.push_back(GridIndexWithVal{0,0});
    gs.push_back(GridIndexWithVal{0,0});
    gs.push_back(GridIndexWithVal{1,0});
    gs.push_back(GridIndexWithVal{0,1});
    gs.push_back(GridIndexWithVal{1,1});
    gs.push_back(GridIndexWithVal{1,2});
    gs.push_back(GridIndexWithVal{2,1});
    gs.push_back(GridIndexWithVal{2,1});
    
    std::set<GridIndexWithVal> set;
    std::set<GridIndexWithVal>::iterator it;
    std::pair<std::set<GridIndexWithVal>::iterator,bool> ret;
    for(GridIndexWithVal& g : gs)
    {   
        ret = set.insert(g);
        if (ret.second == false)
        {
            it=ret.first;
            (*it).val ++;
        }
    }
    
    for(auto g : set)
    {
        g.print();
    }
    
    
    return 0;
}
