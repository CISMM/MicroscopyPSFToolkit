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

#include "xcosm.h"
#include "string.h"
#include "blitz/arrayManip.h"
#include "wu/wuHeader.h"

extern "C" {
    extern float deltar;
    extern float deltaxy;
    extern float deltaxy_nyq;
};

using namespace blitz;

namespace cosm {

template<typename T>
void packHead( RadialPSF<T>* radialPSF,  float* radial, osm_ds& head, int nR, int nZ );

template<typename T>
void CompleteXcosm<T>::rotate(
    bool exact,
    PsfUser* user
){
	std::cout <<"CompleteXcosm::rotate; enter"<<std::endl;
    if ( this->type_ == DIC || this->type_ == DIC_2D )
	{
	    return rotateDIC(user);
    }
    int nR = this->radialPSF_->nR();
    int nZ = this->radialPSF_->nZ();
    int nXY = this->radialPSF_->nXY();
	TinyVector<int, 3> extent(nZ, nXY, nXY);
	
    this->psf_.resize(extent);
	complete_.resize(extent);
    Array<T,2> psf = this->radialPSF_->psf();
    radial_.resize(psf.extent());
    radial_ = cast<float>(psf);
    osm_ds head;
	packHead( this->radialPSF_, radial_.data(),  head, nR, nZ);

    // additional global parameters!    
    deltar = (float) this->radialPSF_->deltaR();
    deltaxy = (float) this->radialPSF_->deltaXY();
    deltaxy_nyq = (float) this->radialPSF_->deltaXYNyq();
	
	double distance = (double)distance_ == 0 || distance_ > 1e15 ? 
        1e20  :							/* single aperture use huge distance */
        distance / (lm_ * magY_ * 1E3);  /* convert to mm in object space */
		
    double fsize = (double)fsize_/(lm_ * magY_ * 1E3); /* convert to mm in object space */
    bool use2Photon =  this->type_ == CONFOCAL_ROTATING_DISK_CIRCULAR_APERTURE && fsize < deltaxy ? true : false;

    fsize = fsize/deltar < 0 ? 1 : fsize/deltar;
	distance = distance/deltar < 0 ? 1 : distance/deltar;	

    switch ( this->type_ )
    {
        case OPTICAL_SECTIONING_WIDEFIELD:
	    rotnone(complete_.data(), &head);
            break;
        case OPTICAL_SECTIONING_2_PHOTON: 
	    rotnone(complete_.data(), &head);
            complete_ = complete_ * complete_;
            break;
        case CONFOCAL_ROTATING_DISK_CIRCULAR_APERTURE:
            if ( use2Photon )
            {
			    rotnone(complete_.data(), &head);
				complete_ = complete_ * complete_;
            }
            else
            {
                rotdiskcirc(complete_.data(),&head,1,distance,fsize,0);   
            }
            break;
        case CONFOCAL_ROTATING_DISK_LINE_APERTURE: 
            rotdiskline(complete_.data(),&head,1,distance,fsize,0); 
            break;
        case DIC:
		case DIC_2D:
            break;
    }
    if ( user != NULL )
    {
            user->update(PsfUser::COMPLETE, nZ, nZ);
    }
    this->psf_ = cast<T>(complete_);
    this->psf_ /= sum(this->psf_);
	std::cout <<"CompleteXcosm::rotate; exit"<<std::endl;
};

template<typename T>
void CompleteXcosm<T>::rotateAndSum(
    bool exact,
    PsfUser* user
){
    std::cout <<"CompleteXcosm::rotateAndSum; enter"<<std::endl;
    if ( this->type_ == DIC )
	{
	    return rotateDIC(user);
    }
    int nR = this->radialPSF_->nR();
    int nZ = this->radialPSF_->nZ();
    int nXY = this->radialPSF_->nXY();
    TinyVector<int, 3> extent(nZ, nXY, nXY);
    this->psf_.resize(extent);
    complete_.resize(extent);
	
    Array<T,2> psf = this->radialPSF_->psf();
    radial_.resize(psf.extent());
    radial_ = cast<float>(psf);
	
    osm_ds head;
	packHead( this->radialPSF_, radial_.data(),  head, nR, nZ);

    // additional global parameters!  
    deltar = (float) this->radialPSF_->deltaR();
    deltaxy = (float) this->radialPSF_->deltaXY();
    deltaxy_nyq = (float) this->radialPSF_->deltaXYNyq();

    double distance = (double)distance_ == 0 || distance_ > 1e15 ? 
        1e20  :							/* single aperture use huge distance */
        distance / (lm_ * magY_ * 1E3);  /* convert to mm in object space */
		
    double fsize = (double)fsize_/(lm_ * magY_ * 1E3); /* convert to mm in object space */
    
    bool use2Photon =  this->type_ == CONFOCAL_ROTATING_DISK_CIRCULAR_APERTURE && fsize < deltaxy ? true : false;

    fsize = fsize/deltar < 0 ? 1 : fsize/deltar;
	distance = distance/deltar < 0 ? 1 : distance/deltar;	

    switch ( this->type_ )
    {
        case OPTICAL_SECTIONING_WIDEFIELD:
            rotsum(complete_.data(), &head);
            break;
        case OPTICAL_SECTIONING_2_PHOTON:
            rotsum(complete_.data(), &head);
            complete_ = complete_ * complete_;
            break;
        case CONFOCAL_ROTATING_DISK_CIRCULAR_APERTURE:
            if ( use2Photon )
            {
			    rotsum(complete_.data(), &head);
				complete_ = complete_ * complete_;
            }
            else
            {
                rdcircsum(complete_.data(),&head,1,distance,fsize,0);   
            }
            break;
        case CONFOCAL_ROTATING_DISK_LINE_APERTURE:
            rdlinesum(complete_.data(),&head,1,distance,fsize,0);
            break;
        case DIC:
		case DIC_2D:
            break;
    }
    if ( user != NULL )
    {
            user->update(PsfUser::COMPLETE, nZ, nZ);
    }

    this->psf_ = cast<T>(complete_);
    this->psf_ /= sum(this->psf_);
    std::cout <<"CompleteXcosm::rotateAndSum; exit"<<std::endl;
};

template<typename T>
void CompleteXcosm<T>::rotateDIC( PsfUser* user )
{
    std::cout <<"CompleteXcosm::rotateDIC; enter"<<std::endl;
    int nR = this->radialPSF_->nR();
    int nZ = this->radialPSF_->nZ();
    int nXY = this->radialPSF_->nXY();
    this->psf_.resize(nZ, nXY, nXY);
	this->psfReal_.resize(nZ, nXY, nXY);
	this->psfImag_.resize(nZ, nXY, nXY);
	std::cout <<"CompleteXcosm::rotateDIC; nR: "<<nR<<", nZ: "<<nZ<<", nXY: "<<nXY << std::endl;
	
	TinyVector<int, 3> extent( nZ, 2*nXY, 2*nXY);
	std::cout <<"CompleteXcosm::rotateDIC; new extent: "<<extent<< std::endl;
    completeRe_.resize(extent);
	completeIm_.resize(extent);

    Array<T,2> psfRe = this->radialPSF_->psfRe();
	Array<T,2> psfIm = this->radialPSF_->psfIm();
    radialRe_.resize(psfRe.extent());
	radialIm_.resize(psfIm.extent());
    radialRe_ = cast<float>(psfRe);
	radialIm_ = cast<float>(psfIm);
    osm_ds headRe;
    osm_ds headIm;
	packHead( this->radialPSF_, radialRe_.data(),  headRe, nR, nZ);
	packHead( this->radialPSF_, radialIm_.data(),  headIm, nR, nZ);
	
    if ( this->type_ == DIC_2D )
    {
        rotinfdicline(completeRe_.data(),completeIm_.data(),&headRe,&headIm,1,shear_,bias_,ampRatio_,0); 
    }
    else
	{
        rotdicline(completeRe_.data(),completeIm_.data(),&headRe,&headIm,1,shear_,bias_,ampRatio_,0); 
    }
    if ( user != NULL )
    {
            user->update(PsfUser::COMPLETE, nZ, nZ);
    }
	wuDataWrite(completeRe_, "dic_comp_re.wu");
	wuDataWrite(completeIm_, "dic_comp_im.wu");
    std::cout <<"CompleteXcosm::rotateDIC; exit"<<std::endl;
}

template<typename T>
void CompleteXcosm<T>::rotateXY( double angle ) 
{
    std::cout <<"CompleteXcosm::rotateXY; enter"<<std::endl;
	int nZ = this->radialPSF_->nZ();
    int nXY = this->radialPSF_->nXY();
	int halfXY = nXY/2;
	int halfZ = nZ/2;
	double cosTheta = cos(angle);
	double sinTheta = sin(angle);
	TinyVector<int,3> shift(halfZ, halfXY, halfXY);
    Array<float,3> shiftedRe = circularShift(completeRe_, shift);
	Array<float,3> shiftedIm = circularShift(completeIm_, shift);
	TinyVector<int,3> extent(nZ, nXY, nXY);
	Array<T,3> newRe(extent);
	Array<T,3> newIm(extent);
	
	for ( int x = 0; x < nXY; x++ )
	{
		for ( int y = 0; y < nXY; y++ )
		{
			double x1 = (x - halfXY) * cosTheta - (y - halfXY) * sinTheta + halfXY;
			double y1 = (x - halfXY) * sinTheta + (y - halfXY) * cosTheta + halfXY;
			int xi = int(x1);
			int yi = int(y1);
			double u = x1 - xi;
			double v = y1 - yi;
			for ( int z = 0; z < nZ; z++ )
			{
				newRe(z, yi, xi) = u*v*shiftedRe(z, yi, xi) + (1-u)*v*shiftedRe(z, yi , xi+1) + 
								   u*(1-v)*shiftedRe( z, yi+1, xi) + (1-u)*(1-v)*shiftedRe(z, yi+1, xi+1);
				newIm(z, yi, xi) = u*v*shiftedIm(z, yi, xi) + (1-u)*v*shiftedIm(z, yi , xi+1) + 
								   u*(1-v)*shiftedIm( z, yi+1, xi) + (1-u)*(1-v)*shiftedIm(z, yi+1, xi+1);
			}
		}
	}
    Array<T,3> psfRe = circularShift(newRe, shift);
	Array<T,3> psfIm = circularShift(newIm, shift);		
	this->psfReal_ = psfRe;
	this->psfReal_ /= sum(this->psfReal_); 
	this->psfImag_ = psfIm;
	this->psfImag_ /= sum(this->psfImag_); 
	this->psf_ = sqrt(psfRe*psfRe + psfIm*psfIm);
	this->psf_ /= sum(this->psf_);
    std::cout <<"CompleteXcosm::rotateXY; exit"<<std::endl;
}

template<typename T>
void packHead( RadialPSF<T>* radialPSF, float* radial, osm_ds& head, int nR, int nZ )
{ 
    memset((void*)&head, 0, 1024);
    head.data = (char*)radial;
    head.mode = WU_FLOAT;
    head.nx = nR; 
    head.ny = nZ;
    head.nz = 1;
    head.user16_footer_size = 0;
    head.user15 = TVAL; 

    head.xstart = radialPSF->oversampling();
    head.ystart = radialPSF->nXY();
    head.zstart = radialPSF->isSymmetric();
    head.ylength = radialPSF->deltaXY();
    head.xlength = radialPSF->deltaR();
}


};
