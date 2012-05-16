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

#ifndef _ESTIMATE_PENALTY_H
#define _ESTIMATE_PENALTY_H

#include <blitz/array.h>

using namespace blitz;

namespace cosm {

template<typename T, int N>
class EstimatePenalty {
   
  public:
    
    EstimatePenalty( T value ) : alpha_(value) {};
    virtual ~EstimatePenalty() {};

    virtual void operator()( Array<T,N>& est ) = 0;   

    T alpha( void ) { return alpha_; };

    void alpha( T value ) { alpha_ = value; }; 

  private: 
    // not allowed
    EstimatePenalty( EstimatePenalty<T,N>& );
    void operator=( EstimatePenalty<T,N>& );

  protected:
     T alpha_;

};

template<typename T, int N>
class IntensityPenalty : public EstimatePenalty<T,N> {
 
  public:

    IntensityPenalty( T value ) : EstimatePenalty<T,N>(value) {};
    ~IntensityPenalty() {};

    void operator()( Array<T,N>& est )   
    {
        est = 1/this->alpha_ *( -1 + sqrt( 1 + 2* this->alpha_ * est));
    };

  private: 
    // not allowed
    IntensityPenalty( IntensityPenalty<T,N>& );
    void operator=( IntensityPenalty<T,N>& );
};

template<typename T, int N>
class RoughnessPenalty : public EstimatePenalty<T,N> {
 
  public:

    RoughnessPenalty( T value ) : EstimatePenalty<T,N>(value) {};
    ~RoughnessPenalty() {};

    void operator()( Array<T,N>& est )   
    {
        Array<T,N> B(est.shape());
        TinyVector<int, N> pos;
        TinyVector<int, N> shift[N];
        TinyVector<int, N> posMin;
        TinyVector<int, N> posMax;
        posMin = est.lbound();;
        posMax = est.ubound();
        T epsilon = T(1e-6); 
        for ( int i = 0; i < est.dimensions(); i++ )
        {
            shift[i] = 0;
            shift[i](i) = 1;
        }

        B = where ( est < epsilon, epsilon, est);
        B = sqrt(B);
        T maxVal = max(B);
        T tol = maxVal/T(20);

        T relax = 1.5;
        T ratdel = T(1);     // NOTE: should be ((deltax+deltay)/(4*deltaz))^2
        T conb1 = T(4) * this->alpha_;
        T conb0 = T(1) +T(2)*conb1*(T(2)+T(ratdel)); 
        T conb = T(1) - relax;
        T conv = T(0.02);
        bool done = false;
        int count = 0;
        while ( !done && (count < 100) )
        {
            T devMax = 0;
            B = where( B < epsilon, epsilon, B );
            ArrayIterator<T,N> iter = B.begin(), end = B.end();
            while ( iter != end ) 
            {
                pos = (*iter);
                T tmpOld = B(pos);
                T tmpDev = -est(pos)/(tmpOld*tmpOld);
                T con2 = relax/(conb0-tmpDev);
                T tmpNew = 0;
                for ( int i = 0; i < est.dimensions(); i++ )
                {
                    TinyVector<int, N> tmp;
                    tmp = min(posMax, max(posMin, pos-shift[i]));
                    tmpNew += B(tmp);
                    tmp = min(posMax, max(posMin, pos+shift[i]));
                    tmpNew += B(tmp);
                } 
                tmpNew = conb*tmpOld + con2*(-2*tmpDev*tmpOld + conb1*tmpNew);
                B(pos) = tmpNew;
                if ( tmpNew > tol ) 
                { 
                    T dev = abs(tmpNew/tmpOld - 1);
                    if ( devMax < dev ) devMax = dev;
                }
                ++iter;
            } 
            if ( devMax <= conv )
            {
                done = true;
            }
            else if ( count == 100 )
            {
                std::cout <<"PDE Roughness penalty has no convergenece"<< std::endl;
            }
            count++;
        }
        est =  pow2(B);
    }

  private: 
    // not allowed
    RoughnessPenalty( RoughnessPenalty<T,N>& );
    void operator=( RoughnessPenalty<T,N>& );
};

};

#endif // _ESTIMATE_PENALTY_H
