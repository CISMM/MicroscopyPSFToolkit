#include "washu.h"

#define MAXPATHLEN 256

extern "C" {

// routine to create volume PSF from XZ cross-section (interp) 
void rotnone(
    float *vol, 
    osm_ds *head
);

// routine to create volume PSF from XZ cross-section (sums nghbr pixels) 
void rotsum(
    float *vol, 
    osm_ds *head
);

// routine to create volume PSF from XZ cross-section for circular confocal 
void rotdiskcirc(
    float *vol, 
    osm_ds *head,
    int bin, 
    float dist, 
    float sz, 
    int tandem
);

// routine to create volume PSF from XZ cross-section for circular confocal 
// by summing neighboring pixels for undersampled case 
void rdcircsum(
    float *vol, 
    osm_ds *head,
    int bin, 
    float dist, 
    float sz, 
    int tandem
);

//routine to create volume PSF from XZ cross-section for slit confocal 
void rotdiskline(
    float *vol, 
    osm_ds *head,
    int bin, 
    float dist, 
    float sz, 
    int tandem
);

//routine to create volume PSF from XY cross-section for slit confocal 
void rotinfdicline(
    float *volRe, 
    float *volIm, 
    osm_ds *headRe,
	osm_ds *headIm,
	int bin, 
    float DeltaX,
    float DeltaPhi,
    float AmpRatio,
    int tandem
);

// routine to create volume PSF from XZ created by gibson_xzslice
void rotdicline(
    float *volRe, 
    float *volIm, 
    osm_ds *headRe,
	osm_ds *headIm,
    int bin, 
    float DeltaX,
    float DeltaPhi,
    float AmpRatio,
    int tandem
);


// routine to create volume PSF from XZ cross-section for slit confocal 
// by summing neighboring pixels for undersampled case 
void rdlinesum(
    float *vol, 
    osm_ds *head,
    int bin, 
    float dist, 
    float sz, 
    int tandem
);



};
