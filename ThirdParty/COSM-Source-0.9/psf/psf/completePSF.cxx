/****************************************************************************
 * Copyright (c) 2004 Einir Valdimarsson and Chrysanthe Preza
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

namespace cosm {

template<typename T>
void CompletePSF<T>::rotate(
    bool exact,
	EvalType type,
    PsfUser* user
){
    Range all = Range::all();
    int nZ = radialPSF_->nZ();
    int nXY = radialPSF_->nXY();
    int nXYHalf = nXY/2;
    int nZHalf = nZ/2;
    Array<T,2> rotatedXY(nXY, nXY);
	switch (type)
	{
        case MAGNITUDE: psf_.resize(nZ, nXY, nXY);
	    case REAL: psfReal_.resize(nZ, nXY, nXY);
	    case IMAGINARY: psfImag_.resize(nZ, nXY, nXY);
	}
    int maxZ = radialPSF_->isSymmetric() ? nZHalf : nZ-1;
    for ( int z = 0; z <= maxZ; z++ ) 
    {   
	for ( int y = 0; y <= nXYHalf; y++ )
	{
	    for ( int x = 0; x <= y; x++ )
	    {
		    T value = (exact ) ? 
			    radialPSF_->exactValue(z,y,x, type) : 
			    radialPSF_->interpolatedValue(z,y,x, type);

		    rotatedXY(y,x) = value;
		    rotatedXY(x,y) = value;

		    if ( x > 0 )
			{
				rotatedXY(y,nXY-x) = value;
  		        rotatedXY(nXY-x,y) = value;
		    }
	 	    if ( y > 0 )
	    	{
		        rotatedXY(x,nXY-y) = value;
		        rotatedXY(nXY-y,x) = value;
		    }
		    if ( y > 0 && x > 0 )
		    {
		        rotatedXY(nXY-y,nXY-x) = value;
		        rotatedXY(nXY-x,nXY-y) = value;
		    }
	    }
	}
	switch (type )
	{
	    case MAGNITUDE:  psf_( z, all, all) = rotatedXY;
		case REAL: psfReal_( z, all, all) = rotatedXY;
		case IMAGINARY: psfImag_( z, all, all) = rotatedXY;
	}
	if ( radialPSF_->isSymmetric() && z > 0 ) 
	{
	    switch (type)
		{
		    case MAGNITUDE: psf_(nZ-z, all, all) = rotatedXY;
	            case REAL: psfReal_(nZ-z, all, all) = rotatedXY;
		    case IMAGINARY: psfImag_(nZ-z, all, all) = 0 - rotatedXY;
		}
	}	
	if ( user != NULL )
	{
	    user->update(PsfUser::COMPLETE, z, maxZ);
	}
    }
//    T maxVal = max(psf_);
//    psf_ /= maxVal;
    switch (type)
	{
	    case MAGNITUDE: psf_ /= sum(psf_);
		default: break;
	}
};

template<typename T>
void CompletePSF<T>::rotateAndSum(
    bool exact,
	EvalType type,
    PsfUser* user
){
    Range all = Range::all();
    int oversampling = exact ? 1 : radialPSF_->oversampling();
    int sq = oversampling * oversampling;
    int nZ = radialPSF_->nZ();
    int nXY = radialPSF_->nXY();
    int nXYHalf = nXY/2;
    int nZHalf = nZ/2;
    Array<T,2> rotatedXY(nXY, nXY);
	switch (type)
	{
        case MAGNITUDE: psf_.resize(nZ, nXY, nXY);
	    case REAL: psfReal_.resize(nZ, nXY, nXY);
	    case IMAGINARY: psfImag_.resize(nZ, nXY, nXY);
	}
    int maxZ = radialPSF_->isSymmetric() ? nZHalf : nZ-1;
    for ( int z = 0; z <= maxZ; z++ ) 
    {   
        for ( int y = 0; y <= nXYHalf; y++ )
        {
            for ( int x = 0; x <= y; x++ )
            {
                T sum = 0;
                for ( int i = 0; i < oversampling; i++ ) 
                {
                    for ( int j = 0; j < oversampling; j++ ) 
                    {
                        sum += exact ?
                            radialPSF_->exactValue(z,y*oversampling+i,x*oversampling+j, type) :
                            radialPSF_->interpolatedValue(z, y*oversampling+i, x*oversampling+j, oversampling, type);
                    }
		        }
		        T value = sum/sq;

		        rotatedXY(y,x) = value;
		        rotatedXY(x,y) = value;
		        if ( x > 0 )
		        {
		            rotatedXY(y,nXY-x) = value;
  		            rotatedXY(nXY-x,y) = value;
		        }
		        if ( y > 0 )
		        {
		            rotatedXY(x,nXY-y) = value;
		            rotatedXY(nXY-y,x) = value;
		        }
		        if ( y > 0 && x > 0 )
		        {
		            rotatedXY(nXY-y,nXY-x) = value;
		            rotatedXY(nXY-x,nXY-y) = value;
		        }
            }
        }
        switch (type )
	    {
	        case MAGNITUDE:  psf_( z, all, all) = rotatedXY;
		    case REAL: psfReal_( z, all, all) = rotatedXY;
		    case IMAGINARY: psfImag_( z, all, all) = rotatedXY;
	    }

	    if ( radialPSF_->isSymmetric() && z > 0 ) 
	    {
	        switch (type)
		    {
		        case MAGNITUDE: psf_(nZ-z, all, all) = rotatedXY;
	            case REAL: psfReal_(nZ-z, all, all) = rotatedXY;
			    case IMAGINARY: psfImag_(nZ-z, all, all) = 0 - rotatedXY;
		    }
	    }
	    if ( user != NULL )
	    {
	        user->update(PsfUser::COMPLETE, z, maxZ);
	    }
    }
//    T maxVal = (max)(psf_);
//    psf_ /= maxVal;
      switch (type)
      {
	    case MAGNITUDE: psf_ /= sum(psf_);
		default: break;
      }
};

template<typename T>
void CompletePSF<T>::sumY(
    void
) {
    Range all = Range::all();
    int nZ = radialPSF_->nZ();
    int nXY = radialPSF_->nXY();
    firstIndex i;
    secondIndex j;
    thirdIndex k;
    Array<T,2> psfReduced;
    psfReduced.resize(nZ,nXY);
    psfReduced  = sum(psf_(i,k,j),k); 
    psf_.resize(nZ, 1, nXY);
    psf_(all,0,all) = psfReduced;
}


};
