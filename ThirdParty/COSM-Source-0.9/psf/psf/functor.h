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
                                                                                
#ifndef _FUNCTOR_H
#define _FUNCTOR_H

namespace cosm {

template<typename T>
class Functor {

  public:

    Functor() {};
    virtual ~Functor() {};

    virtual T operator()( T /*x*/ ) { return 0; };
    virtual T operator()( T /*x*/, T /*y*/ ) { return 0; };
    virtual T operator()( T /*x*/, T /*y*/, T /*z*/ ) { return 0; };
  
  protected:

    // not allowed
    Functor( Functor<T>& );
    Functor& operator=( Functor<T>& );

};
};

#endif // _FUNCTOR_H
