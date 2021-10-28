#pragma once
#include "Eigen/Dense"
#include <inttypes.h>
#include <cmath>
#include <vector>
#include "mshw_quatmath.h"

namespace ocgrid_constants
{
    static const float eps = 1e-3;
};

struct GridIndex
{
    int x;
    int y;
    
    void print() const
    {
        std::cout << "x: " << x << " y: " << y << std::endl; 
    };
};

struct NeoSector
{
    Eigen::Vector2f sector_origin;
    Eigen::Vector2f sector_dir;
    float sector_angle;
    float sector_range;
    float mes_range;
    float free_addition;
    float busy_addition;
};

enum MapLayer
{
    COMMON,
    HUMAN,
    PILLAR,
    VEHICLE,
    ANIMAL,
    OTHER,
    UNKNOWN,
    MAP_LAYER_END
};

class OcGrid
{
public:
    OcGrid(float min_size, float resolution, float min_free_val, float max_free_val, float min_busy_val, float max_busy_val)
    {
        //_human_ocmap; _pillar_ocmap; _vehicle_ocmap; _animal_ocmap; _other_ocmap;
        _resolution = resolution;
        _dimention = (int)std::ceil(min_size/resolution);
        _size = _dimention * _resolution;
        _min_free_val = min_free_val;
        _max_busy_val = max_busy_val;
        
        for (int i = COMMON; i != MAP_LAYER_END; i++)
        {
            Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic> map;
            map.setConstant(_dimention, _dimention, 0);
            _ocmap_layers.push_back(map);
        }
    };
            
    void refresh_grid_with_lim(const std::vector<Eigen::Vector2f>& p1s, const std::vector<Eigen::Vector2f>& p2s, float lim, float free_addition, float busy_addition, MapLayer map_layer)
    {
        for (int i = 0; i < p1s.size(); i ++)
        {
            Eigen::Vector2f p1 = p1s.at(i);
            Eigen::Vector2f p2 = p2s.at(i);

            std::vector<GridIndex> on_ray_opened = onray_indexes_with_lim(p1, p2, lim);
            for (const GridIndex& index : on_ray_opened)
            {
                std::cout << "x: " << index.x << " y: " << index.y << std::endl;
                add_value_on_index_with_lim_check(index, free_addition, map_layer);
            }
        }
                
        for (int i = 0; i < p1s.size(); i ++)
        {
            Eigen::Vector2f p1 = p1s.at(i);
            Eigen::Vector2f p2 = p2s.at(i);

            std::vector<GridIndex> end_ray_closed =  endray_indexes_with_lim(p1, p2, lim);
            for (const GridIndex& index : end_ray_closed)
            {
                add_value_on_index_with_lim_check(index, busy_addition, map_layer);
            }
        }
    };
    
    void refresh_grid_with_sector(const std::vector<Eigen::Vector2f>& p1s, const std::vector<Eigen::Vector2f>& p2s,
                                  const Eigen::Vector3f& r, const Eigen::Vector4f& q,
                                  float r_lim, float phi_lim, 
                                  float free_addition, float busy_addition,
                                  MapLayer map_layer)
    {
        Eigen::Vector2f sector_origin(r(0), r(1));
        Eigen::Vector3f dir3 = quatRotate(q, Eigen::Vector3f(1,0,0));
        Eigen::Vector2f sector_dir(dir3(0), dir3(1));
        std::vector<GridIndex> sector_opened =  sector_internal_indexes(sector_origin, sector_dir, r_lim, phi_lim);
        
        //std::cout << "sector_origin: " << sector_origin.transpose() << std::endl;
        //std::cout << "sector_dir: " << sector_dir.transpose() << std::endl;
        //std::cout << "so_size: " << sector_opened.size() << std::endl;
        //std::cout << "p1s_size: " << p1s.size() << std::endl;
        
        for (const GridIndex& index : sector_opened)
        {
            add_value_on_index_with_lim_check(index, free_addition, map_layer);
            //std::cout << "x : " << index.x << std::endl;
            //std::cout << "y : " << index.y << std::endl;
            //index.print();
        }
                
        for (int i = 0; i < p1s.size(); i ++)
        {
            Eigen::Vector2f p1 = p1s.at(i);
            Eigen::Vector2f p2 = p2s.at(i);

            std::vector<GridIndex> end_ray_closed =  endray_indexes_with_lim(p1, p2, r_lim);
            for (const GridIndex& index : end_ray_closed)
            {
                add_value_on_index_with_lim_check(index, busy_addition, map_layer);
                //std::cout << "a : " << index.x << std::endl;
                //std::cout << "b : " << index.y << std::endl;
            }
        }
    };
    
    void refresh_grid_with_neo_sectors(const std::vector<NeoSector>& neo_sectors, MapLayer map_layer)
    {
        
        
        /* opened */
        /*
        for (auto& s : neo_sectors)
        {
            std::vector<GridIndex> opened = 
                sector_internal_indexes(s.sector_origin, s.sector_dir, s.sector_range, s.sector_angle);
            for (auto& o : opened)
            {
                add_value_on_index_with_lim_check(o, s.free_addition, map_layer);
            }
        }
        */
        
        /* closed */
        for (auto& s : neo_sectors)
        {
            std::vector<GridIndex> closed;
            // check, if sector so small we use ray update procedure to avoid empty cell sector
            if (s.mes_range * s.sector_angle > 2*_resolution && s.mes_range > 2*_resolution)
            {
                closed = sector_border_indexes(s.sector_origin, s.sector_dir, s.mes_range, s.sector_angle);
            }
            else
            {
                closed = endray_indexes_with_lim(
                    s.sector_origin, s.sector_origin + s.sector_dir * s.mes_range, s.sector_range);
            }
            for (auto& c : closed)
            {
                add_value_on_index_with_lim_check(c, s.busy_addition, map_layer);
            }
        }
    };
    
    const Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic>& get_map(MapLayer map_layer)
    {
        return _ocmap_layers.at(map_layer);
    };
    
    std::vector<Eigen::Vector2f> get_obstcl_points(float busy_lim, MapLayer map_layer)
    {
        std::vector<Eigen::Vector2f> res;
        std::vector<GridIndex> inds = find_cells_more_than(busy_lim, map_layer);
        for (GridIndex index : inds)
        {
            res.push_back(center_of_cell(index));
        }
        return res;
    };
    
    ~OcGrid()
    {
    };

    void clear()
    {
        for (int i = COMMON; i != MAP_LAYER_END; i++)
        {
            _ocmap_layers.at(i).setConstant(_dimention, _dimention, 0);
        }
    }
    
    void exponentially_foget(float rate)
    {
        
        for (int i = COMMON; i != MAP_LAYER_END; i++)
        {
            Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic>& map = _ocmap_layers.at(i);
            map = (map.cast<float>() * rate).cast<int8_t>();
        }
    }
    
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
    
    void print_map(MapLayer map_layer)
    {
        std::cout << std::endl << get_map(map_layer).transpose().cast<int>() << std::endl;
    }
    
private:
    bool check_index_in(const GridIndex& index)
    {
        if (index.x >= _dimention || index.y >= _dimention || index.x < 0 || index.y < 0)
        {
            return false;
        }
        return true;
    };
    
    bool check_point_in(const Eigen::Vector2f& p)
    {
        float grid_lim = _size / 2.f - ocgrid_constants::eps;
        if (p.x() > grid_lim || p.x() < -grid_lim || p.y() > grid_lim || p.y() < -grid_lim)
        {
            return false;
        }
        return true;
    };
    
    std::vector<Eigen::Vector2f> x_area(const Eigen::Vector2f& p, float range)
    {
        
        Eigen::Vector2f p_mp = p + Eigen::Vector2f(-range, range);
        Eigen::Vector2f p_pp = p + Eigen::Vector2f(range, range);
        Eigen::Vector2f p_mm = p + Eigen::Vector2f(-range, -range);
        Eigen::Vector2f p_pm = p + Eigen::Vector2f(range, -range);
        
        std::vector<Eigen::Vector2f> points;
        points.push_back(p_mp);
        points.push_back(p_pp);
        points.push_back(p_mm);
        points.push_back(p_pm);
        return points;
    };
    
    std::vector<Eigen::Vector2f> h_area(const Eigen::Vector2f& p, float range)
    {
        Eigen::Vector2f p_p = p + Eigen::Vector2f(range, 0);
        Eigen::Vector2f p_m = p + Eigen::Vector2f(-range, 0);

        std::vector<Eigen::Vector2f> points;
        points.push_back(p_p);
        points.push_back(p_m);
        return points;
    };
    
    std::vector<Eigen::Vector2f> v_area(const Eigen::Vector2f& p, float range)
    {
        Eigen::Vector2f p_p = p + Eigen::Vector2f(0, range);
        Eigen::Vector2f p_m = p + Eigen::Vector2f(0, -range);

        std::vector<Eigen::Vector2f> points;
        points.push_back(p_p);
        points.push_back(p_m);
        return points;
    };
    
    std::vector<float> linespace(float min, float step, float max)
    {
        std::vector<float> res;
        float value = min;
        while (value <= max)
        {
            res.push_back(value);
            value += step;
        }
        return res;
    }
    
    GridIndex index_of_point(const Eigen::Vector2f& p)
    {
        GridIndex index;
        index.x = std::floor(p.x() / _resolution + (float)_dimention / 2.f);
        index.y = std::floor(p.y() / _resolution + (float)_dimention / 2.f);
        return index;
    };
    
    void set_value_on_index(const GridIndex& index, int8_t value, MapLayer map_layer)
    {
        _ocmap_layers.at(map_layer)(index.x,index.y) = value;
    }
    
    void add_value_on_index_with_lim_check(const GridIndex& index, int8_t value, MapLayer map_layer)
    {   
        Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic>& map = _ocmap_layers.at(map_layer);
        map(index.x,index.y) += value;
        if (map(index.x,index.y) > _max_busy_val)   
        {
            map(index.x,index.y) = _max_busy_val;
        }
        if (map(index.x,index.y) < _min_free_val)
        {
            map(index.x,index.y) = _min_free_val;
        }
    }
        
    std::vector<GridIndex> x_area_indexes_of_point(const Eigen::Vector2f& p)
    {
        std::vector<GridIndex> indexes;
        if (!check_point_in(p))
        {
            return indexes;
        }
        
        std::vector<Eigen::Vector2f> x_area_poitns = x_area(p, ocgrid_constants::eps);
        for (const Eigen::Vector2f& x_point : x_area_poitns)
        {
            GridIndex index = index_of_point(x_point);
            if (check_index_in(index))
            {
                indexes.push_back(index);
            }
        }
        return indexes;
    };
    
    std::vector<GridIndex> endray_indexes_with_lim(const Eigen::Vector2f& p1, const Eigen::Vector2f& p2, float lim)
    {
        float eps = ocgrid_constants::eps;
        std::vector<GridIndex> indexes;
        Eigen::Vector2f dp = p2 - p1;
        float ndp = dp.norm();
        if (ndp < eps)
        {
            return indexes;
        }
            
        if (ndp > lim)
        {
            return indexes;
        }
        
        return x_area_indexes_of_point(p2);
    }
    
    std::vector<GridIndex> onray_indexes_with_lim(const Eigen::Vector2f& p1, const Eigen::Vector2f& p2, float lim)
    {
        float eps = ocgrid_constants::eps;
        std::vector<GridIndex> indexes;
        Eigen::Vector2f dp = p2 - p1;
        float ndp = dp.norm();
        if (ndp < eps)
        {
            return indexes;
        }
        Eigen::Vector2f dir = dp / ndp;
        
        if (ndp > lim)
        {
            ndp = lim;
        }
        
        Eigen::Vector2f p2_cut = p1 + dir*ndp;
        return onray_indexes(p1, p2_cut);
    };
    
    std::vector<GridIndex> sector_internal_indexes(const Eigen::Vector2f& sector_origin,
                                                   const Eigen::Vector2f& sector_dir,
                                                   float r_lim, float phi_lim)
    {
        float eps = ocgrid_constants::eps;
        std::vector<GridIndex> indexes;
        
        Eigen::Vector4f q_half_phi(cos(phi_lim/4), 0, 0, sin(phi_lim/4));
        Eigen::Vector2f a = sector_dir * r_lim;
        float a_phi = atan2(a(1), a(0));
        Eigen::Vector3f a3(a(0), a(1), 0);
        Eigen::Vector3f b3 = quatRotate(q_half_phi, a3);
        Eigen::Vector3f c3 = quatRotate(quatInverse(q_half_phi), a3);
        Eigen::Vector2f b(b3(0), b3(1));
        Eigen::Vector2f c(c3(0), c3(1));
        Eigen::Vector2f a_end = sector_origin + a;
        Eigen::Vector2f b_end = sector_origin + b;
        Eigen::Vector2f c_end = sector_origin + c;
        
        float xmin = std::min( std::min(sector_origin(0), a_end(0)),  std::min(b_end(0), c_end(0)) );
        float xmax = std::max( std::max(sector_origin(0), a_end(0)),  std::max(b_end(0), c_end(0)) );
        float ymin = std::min( std::min(sector_origin(1), a_end(1)),  std::min(b_end(1), c_end(1)) );
        float ymax = std::max( std::max(sector_origin(1), a_end(1)),  std::max(b_end(1), c_end(1)) );
        
        xmin = std::max(-_size/2.f + eps, xmin);
        xmax = std::min(_size/2.f - eps, xmax);
        ymin = std::max(-_size/2.f + eps, ymin);
        ymax = std::min(_size/2.f - eps, ymax);
        
        int xmin_ind = static_cast<int>((xmin + _size/2) / _resolution);
        int xmax_ind = static_cast<int>((xmax + _size/2) / _resolution);
        int ymin_ind = static_cast<int>((ymin + _size/2) / _resolution);
        int ymax_ind = static_cast<int>((ymax + _size/2) / _resolution);

        for (int x_ind = xmin_ind; x_ind <= xmax_ind; x_ind++)
        {
            for (int y_ind = ymin_ind; y_ind <= ymax_ind; y_ind++)
            {
                Eigen::Vector2f cell_center = center_of_cell(GridIndex {x_ind, y_ind});
                Eigen::Vector2f ray = cell_center - sector_origin;
                float dist = ray.norm();
                float phi = atan2(ray(1), ray(0));
                float abs_rot = std::abs(shortestRotation(a_phi, phi));
                if (dist < r_lim && abs_rot < phi_lim / 2)
                {
                    GridIndex index{x_ind, y_ind};
                    if (check_index_in(index))
                    {
                        indexes.push_back(index);
                    }
                }
            }
        }
        
        return indexes;
        //
    };
    
    std::vector<GridIndex> sector_border_indexes(const Eigen::Vector2f& sector_origin,
                                                 const Eigen::Vector2f& sector_dir,
                                                 float r_lim, float phi_lim)
    {
        float eps = ocgrid_constants::eps;
        std::vector<GridIndex> indexes;
        
        Eigen::Vector4f q_half_phi(cos(phi_lim/4), 0, 0, sin(phi_lim/4));
        Eigen::Vector2f a = sector_dir * r_lim;
        float a_phi = atan2(a(1), a(0));
        Eigen::Vector3f a3(a(0), a(1), 0);
        Eigen::Vector3f b3 = quatRotate(q_half_phi, a3);
        Eigen::Vector3f c3 = quatRotate(quatInverse(q_half_phi), a3);
        Eigen::Vector2f b(b3(0), b3(1));
        Eigen::Vector2f c(c3(0), c3(1));
        Eigen::Vector2f a_end = sector_origin + a;
        Eigen::Vector2f b_end = sector_origin + b;
        Eigen::Vector2f c_end = sector_origin + c;
        
        float xmin = std::min( std::min(sector_origin(0), a_end(0)),  std::min(b_end(0), c_end(0)) );
        float xmax = std::max( std::max(sector_origin(0), a_end(0)),  std::max(b_end(0), c_end(0)) );
        float ymin = std::min( std::min(sector_origin(1), a_end(1)),  std::min(b_end(1), c_end(1)) );
        float ymax = std::max( std::max(sector_origin(1), a_end(1)),  std::max(b_end(1), c_end(1)) );
        
        xmin = std::max(-_size/2.f + eps, xmin);
        xmax = std::min(_size/2.f - eps, xmax);
        ymin = std::max(-_size/2.f + eps, ymin);
        ymax = std::min(_size/2.f - eps, ymax);
        
        int xmin_ind = static_cast<int>((xmin + _size/2) / _resolution);
        int xmax_ind = static_cast<int>((xmax + _size/2) / _resolution);
        int ymin_ind = static_cast<int>((ymin + _size/2) / _resolution);
        int ymax_ind = static_cast<int>((ymax + _size/2) / _resolution);

        for (int x_ind = xmin_ind; x_ind <= xmax_ind; x_ind++)
        {
            for (int y_ind = ymin_ind; y_ind <= ymax_ind; y_ind++)
            {
                Eigen::Vector2f cell_center = center_of_cell(GridIndex {x_ind, y_ind});
                Eigen::Vector2f ray = cell_center - sector_origin;
                float dist = ray.norm();
                float phi = atan2(ray(1), ray(0));
                float abs_rot = std::abs(shortestRotation(a_phi, phi));
                float d_res = sqrt(2.f)*_resolution / 2.f;
                if ((dist > (r_lim - d_res)) && (dist < (r_lim + d_res)) && (abs_rot < phi_lim / 2.f))
                {
                    GridIndex index{x_ind, y_ind};
                    if (check_index_in(index))
                    {
                        indexes.push_back(index);
                    }
                }
            }
        }
        
        return indexes;
        //
    };
    
    std::vector<GridIndex> onray_indexes(const Eigen::Vector2f& p1, const Eigen::Vector2f& p2)
    { 
        float eps = ocgrid_constants::eps;
        std::vector<GridIndex> indexes;
        
        Eigen::Vector2f dp = p2 - p1;
        float ndp = dp.norm();
        if (ndp < eps)
        {
            return indexes;
        }
        Eigen::Vector2f dir = dp / ndp;
        
        /* find searching borders */
        float xmin = std::min(p1.x(), p2.x());
        float xmax = std::max(p1.x(), p2.x());
        float ymin = std::min(p1.y(), p2.y());
        float ymax = std::max(p1.y(), p2.y());

        xmin = std::max(-_size/2.f+ eps, xmin);
        xmax = std::min(_size/2.f- eps, xmax);
        ymin = std::max(-_size/2.f + eps, ymin);
        ymax = std::min(_size/2.f - eps, ymax);

        float xmin_d = std::floor((xmin + _size/2) / _resolution) * _resolution - _size/2;
        float xmax_d = std::ceil((xmax + _size/2) / _resolution) * _resolution - _size/2;
        float ymin_d = std::floor((ymin + _size/2) / _resolution) * _resolution - _size/2;
        float ymax_d = std::ceil((ymax + _size/2) / _resolution) * _resolution - _size/2;
        
        /* find intrscts */
        std::vector<float> xrays = linespace(xmin_d + _resolution, _resolution, xmax_d - _resolution);
        std::vector<float> yrays = linespace(ymin_d + _resolution, _resolution, ymax_d - _resolution);
                
        std::vector<Eigen::Vector2f> x_inters_points;
        if (std::abs(dir.x()) > eps)
        {
            for (float x : xrays)
            {
                float t = (x - p1.x()) / dir.x();
                if (t > 0 && t < ndp)
                {
                    float y = p1.y() + dir.y()*t;
                    x_inters_points.push_back(Eigen::Vector2f(x,y));
                }
            }
        }

        std::vector<Eigen::Vector2f> y_inters_points;
        if (std::abs(dir.y()) > eps)
        {
            for (float y : yrays)
            {
                float t = (y - p1.y()) / dir.y();
                if (t > 0 && t < ndp)
                {
                    float x = p1.x() + dir.x()*t;
                    y_inters_points.push_back(Eigen::Vector2f(x,y));
                }
            }
        }
        
        /* find cells */
        for (const Eigen::Vector2f& p : x_inters_points)
        {
            std::vector<Eigen::Vector2f> h_area_poitns = h_area(p, ocgrid_constants::eps);
            for (const Eigen::Vector2f& h_point : h_area_poitns)
            {
                GridIndex index = index_of_point(h_point);
                if (check_index_in(index))
                {
                    indexes.push_back(index);
                }
            }
        }
        
        for (const Eigen::Vector2f& p : y_inters_points)
        {
            std::vector<Eigen::Vector2f> v_area_poitns = v_area(p, ocgrid_constants::eps);
            for (const Eigen::Vector2f& v_point : v_area_poitns)
            {
                GridIndex index = index_of_point(v_point);
                if (check_index_in(index))
                {
                    indexes.push_back(index);
                }
            }
        }
        
        return indexes;
    };
        
    std::vector<GridIndex> find_cells_more_than(int8_t val, MapLayer map_layer)
    {
        std::vector<GridIndex> res;
        Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic>& map = _ocmap_layers.at(map_layer);
        
        for (int x = 0; x < _dimention; x++)
        {
            for (int y = 0; y < _dimention; y++)
            {
                if (map(x,y) > val)
                {
                    GridIndex gi;
                    gi.x = x;
                    gi.y = y;
                    res.push_back(gi);
                }
            }
        }
        
        return res;
    };
    
    Eigen::Vector2f center_of_cell(GridIndex index)
    {
        float center_x = -(_resolution * _dimention)/2 + (_resolution * index.x) + _resolution/2;
        float center_y = -(_resolution * _dimention)/2 + (_resolution * index.y) + _resolution/2;
        return Eigen::Vector2f(center_x, center_y);
    }

private:
    std::vector<Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic>> _ocmap_layers;
    float _size;
    float _resolution;
    int _dimention;
    
    float _min_free_val;
    float _max_busy_val;
};
