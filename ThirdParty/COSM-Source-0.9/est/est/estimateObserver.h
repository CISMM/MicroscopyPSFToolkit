/****************************************************************************
 * Copyright (c) 2007 Einir Valdimarsson and Chrysanthe Preza
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 ****************************************************************************/

#ifndef _ESTIMATE_OBSERVER_H
#define _ESTIMATE_OBSERVER_H

namespace cosm {

template<typename T, int N>
class EstimateIterative;

template<typename T, int N>
class EstimateObserver {
   
  public:
    
    EstimateObserver( unsigned int value ) : period_(value), count_(value) {};
    virtual ~EstimateObserver() {};
    virtual void notify( EstimateIterative<T,N>& estimate ) 
    {
       if ( count_ == 0 ) 
       {
           update(estimate);
           count_ = period_;
       }
       count_--; 
    };

  protected: 
    virtual void update( EstimateIterative<T,N>& estimate ) = 0;

  private: 
    // not allowed
    EstimateObserver( EstimateObserver<T,N>& );
    void operator=( EstimateObserver<T,N>& );

  protected:
     unsigned int period_;
     unsigned int count_; 

};

};

#endif // _ESTIMATE_OBSERVER_H
