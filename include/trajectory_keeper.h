#pragma once
#include "Eigen/Dense"
#include <inttypes.h>
#include <cmath>
#include <vector>
#include <stdio.h>

struct TrajectoryKeeperPoint
{
    double t;
    Eigen::Vector3f r;
    Eigen::Vector4f q;
    void print()
    {
        std::cout << "t: " << t << " r: " << r.transpose() << " q: " << q.transpose() << std::endl;
    };
};

class TrajectoryKeeper
{
public:
    TrajectoryKeeper()
    {
    };
    
    void reinit(int buffer_size)
    {
        _buffer_size = buffer_size;
        _buffer = std::vector<TrajectoryKeeperPoint>();
        _buffer.reserve(_buffer_size + 1);
        _is_full = false;
        _skip_counter = 0;
    };
    
    bool addPoint(double t, const Eigen::Vector3f& r,const Eigen::Vector4f& q)
    {   

        if(_buffer.size() > 0)
        {
            double eps = 1e-6;
            if(t < _buffer.back().t + eps)
            {
                //std::cout << 
                //"ocgrid trajectory keeper error: t_storage < _buffer.back().t + eps, skip traj point"
                //<< std::endl;
                
                //std::cout << 
                //"t0: " << _buffer.back().t << " t1: " << t
                //<< std::endl;
                
                _skip_counter ++;
                
                if (_skip_counter > 5)
                {
                    //std::cout << 
                    //"ocgrid trajectory keeper error: _skip_counter > 5, reinit"
                    //<< std::endl;
                    reinit(_buffer_size);
                }
                return _is_full;
            }
        }
        
        TrajectoryKeeperPoint p = {t, r, q};
        _buffer.push_back(p);
        _skip_counter = 0;
        //_buffer.back().print();
        if(_buffer.size() > _buffer_size)
        {
            _buffer.erase(_buffer.begin());
            _is_full = true;
            
            /*
            for(int i = 0; i < _buffer_size; i++)
            {
                _buffer.at(i).print();
            }
            */
            
        }

        return _is_full;
    };
    
    bool getApproximation(double t, TrajectoryKeeperPoint& p_out)
    {
        bool found = false;
        
        if (!_is_full)
        {
            return found;
        }
        
        if (t < _buffer.front().t)
        {
            return found;
        }
        
        if (t > _buffer.back().t)
        {
            if (t > _buffer.back().t + 0.5)
            {
                //std::cout << "ocgrid trajectory keeper error: t_approx > _buffer.back().t + 0.5s" << std::endl;
                return found;
            }
            else
            {
                p_out = _buffer.back();
                found = true;
                return found;
            }
        }
        
        for(int i = 1; i < _buffer_size; i++)
        {
            TrajectoryKeeperPoint p2 = _buffer.at(i);
            if (t < p2.t)
            {
                TrajectoryKeeperPoint p1 = _buffer.at(i-1);
                double dt = p2.t - p1.t;
                Eigen::Vector3f dr = p2.r - p1.r;
                Eigen::Vector4f dq = p2.q - p1.q;
                double a = t - p1.t;
                
                p_out.t = t;
                p_out.r = p1.r + dr/dt * a;
                p_out.q = p1.q + dq/dt * a;
                p_out.q = p_out.q / p_out.q.norm();
                found = true;
                break;
            }
        }
        return found;
    };
    
bool is_full()
{
    return _is_full;
};
    
private:
    int _buffer_size;
    std::vector<TrajectoryKeeperPoint> _buffer;
    bool _is_full;
    int _skip_counter;
};
