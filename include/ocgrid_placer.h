#pragma once
#include "Eigen/Dense"
#include "mshw_ocgrid.h"
#include <cmath>

class OcGridPlacer
{
public:
    
    OcGridPlacer()
    {
    };
    
    void reinit(const OcGrid& ocGrid, const Eigen::Vector2f& origin)
    {
        _size = ocGrid.get_size();
        _resolution = ocGrid.get_resolution();
        _dimention = ocGrid.get_dimention();
        _origin = origin;
        
    };
    
    bool replaceMapProcedure(const Eigen::Vector2f r_local, OcGrid& ocGrid)
    {
        bool map_replaced = false;
        
        Eigen::Vector2f r_map = localToOcGridFrame(r_local);
        GridIndex ind = ocgrid_index_of_point(r_map, _dimention, _resolution);
        if (fabs(r_map.x()) > _size / 4 || fabs(r_map.y()) > _size / 4)
        {
            map_replaced = true;
            GridIndex ind_center = ocgrid_index_of_point(Eigen::Vector2f(0,0), _dimention, _resolution);
            GridIndex ind_current = ocgrid_index_of_point(r_map, _dimention, _resolution);
            int x_cell_shift = ind_current.x - ind_center.x;
            int y_cell_shift = ind_current.y - ind_center.y;
            if (abs(x_cell_shift) > _dimention-1 || abs(y_cell_shift) > _dimention-1)
            {
                while(true) // for stop send ocgrid data
                {
                    std::cout << "WrOcgrid: ocgrid placer: ERROR: CAN'T MAKE REPLACE MAP PROCEDURE: WRONG MAP SHIFT" << std::endl;
                    usleep(1000000);
                }
            }
            
            Eigen::Vector2f dr_shift(x_cell_shift*_resolution, y_cell_shift*_resolution);
            _origin = _origin + dr_shift;
            ocGrid.shift_map(-x_cell_shift, -y_cell_shift);
        }
        return map_replaced;
    };
    
    Eigen::Vector2f localToOcGridFrame(const Eigen::Vector2f r)
    {
        return r - _origin;
    };
    
    Eigen::Vector2f ocGridToLocalFrame(const Eigen::Vector2f r)
    {
        return r + _origin;
    };
        
    Eigen::Vector2f get_origin()
    {
        return _origin;
    };
    
    float get_size()
    {
        return _size;
    };
    
    float get_resolution()
    {
        return _resolution;
    };
    
    int get_dimention()
    {
        return _dimention;
    };
        
private:
    Eigen::Vector2f _origin;
    float _size;
    float _resolution;
    int _dimention;
};
