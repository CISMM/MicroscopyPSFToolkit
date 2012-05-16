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
 
#ifndef _PSF_FUNCTOR_H
#define _PSF_FUNCTOR_H

#include "psf/functorComplex.h"

namespace cosm {

template<typename T>
class PsfFunctor : public FunctorComplex<T> {

  public:

    PsfFunctor() {};
    ~PsfFunctor() {};

    virtual bool isSymmetric() { return false; };

  protected:

    // not allowed
    PsfFunctor( PsfFunctor<T>& );
    PsfFunctor& operator=( PsfFunctor<T>& );

};
};

#endif // _PSF_FUNCTOR_H
