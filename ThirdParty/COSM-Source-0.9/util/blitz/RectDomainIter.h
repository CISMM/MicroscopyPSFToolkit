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

#ifndef _RECT_DOMAIN_ITER_H
#define _RECT_DOMAIN_ITER_H

#include <blitz/array/domain.h>
#include <blitz/tinyvec.h>

BZ_NAMESPACE(blitz)

template<int N>
class RectDomainIter {

  public:

    RectDomainIter( 
	RectDomain<N>& rectDomain 
    ) : count_(1), rectDomain_(rectDomain), 
	current_(rectDomain.lbound()), end_(0x7fffffff)
    { 
	for ( int i = 0; i < N; i++ ) 
	{
	   count_ *=  (rectDomain_.ubound(i) - rectDomain_.lbound(i) + 1); 	
	   if ( count_ < 0 ) {
		count_ = 0;
	   }
	}
    };

    const RectDomain<N>& domain() const
    {   
	return rectDomain_; 
    };

    const TinyVector<int,N>& current() const
    {   
	return current_; 
    };

    const TinyVector<int,N>& begin() const
    {   
	return (count_ == 0 ) ? end_ : rectDomain_.lbound(); 
    };

    const TinyVector<int,N>& end() const
    {   
	return end_;
    };

    const TinyVector<int,N>& next()
    {   
	if ( count_ != 0 && adjust(N-1) == true ) 
	{ 
	    return current_;
	}
	return end_;
    };

    int count() const
    {
	return count_;
    };

  protected:

    bool adjust(int index) 
    {  
	current_(index) = current_(index) + 1; 
	if ( index == 0 && current_(0) > rectDomain_.ubound(0) ) 
	{
	    current_(0) = rectDomain_.ubound(0);
	    return false;
	}
	if ( index != 0 && current_(index) > rectDomain_.ubound(index) )
        {
	    current_(index) = rectDomain_.lbound(index);
	    return adjust(index-1);
        }
        return true;
    };

  private:

    int count_;
    RectDomain<N> rectDomain_;
    TinyVector<int, N> current_;
    TinyVector<int,N> end_; 

};

template<class T, int N>
bool equal(
    const TinyVector<T,N>& u,
    const TinyVector<T,N>& v
) {
    for ( int i = 0; i < N; i++ ) 
    {
	if ( u(i) != v(i) ) 
	{
	     return false;
	}
    }
    return true;
}

template<class T, int N>
bool larger(
    const TinyVector<T,N>& u,
    const TinyVector<T,N>& v
) {
    for ( int i = 0; i < N; i++ ) 
    {
	if ( u(i) < v(i) ) 
	{
	     return false;
	}
    }
    return true;
}

template<class T, int N>
TinyVector<T,N> (max)(
    const TinyVector<T,N>& u,
    const TinyVector<T,N>& v
) {
    TinyVector<T,N> w;
    for ( int i = 0; i < N; i++ ) 
    {
	w(i) = 	u(i) > v(i) ? u(i) : v(i); 
    }
    return w;
}

template<class T, int N>
TinyVector<T,N> (min)(
    const TinyVector<T,N>& u,
    const TinyVector<T,N>& v
) {
    TinyVector<T,N> w;
    for ( int i = 0; i < N; i++ ) 
    {
	w(i) = 	u(i) < v(i) ? u(i) : v(i); 
    }
    return w;
}

BZ_NAMESPACE_END

#endif // _RECT_DOMAIN_ITER_H
