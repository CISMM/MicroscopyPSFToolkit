/***************************************************************************
  COPYRIGHT 1996 BIOMEDICAL COMPUTER LABORATORY, WASHINGTON UNIVERSITY, MO
****************************************************************************/
/************************************************************************ 

   WRITE_DATASET.C	Procedure for writing a dataset with a 
			WASHU header.

   Author:
	Keith Doolittle

************************************************************************/
#include <stdio.h>
#include <sys/types.h>
#ifndef WIN32
#include <sys/param.h>
#endif
#include "washu.h"

int write_dataset (osm_ds *ds, char *filename)
{
FILE *fp;
int  data_size;

if((ds==(osm_ds*)NULL)||(ds->user15 != TVAL)) {
   fprintf(stderr,"ERROR: Invalid header in write_dataset()\n");
   return(1);
   }

if((ds->nx < 0)||(ds->nx > WU_MAXDIM)||(ds->ny < 0)||(ds->ny > WU_MAXDIM)||
   (ds->nz < 0)||(ds->nz > WU_MAXDIM)) {
   fprintf(stderr,"ERROR: Dimensions in header invalid in write_dataset()\n");
   fprintf(stderr,"       Image dimensions must be in range [0..%d]\n",
					WU_MAXDIM);
   return(1);
   }

if((fp=fopen(filename,"wb"))==(FILE*)NULL) {
   fprintf(stderr,"ERROR: Can't open `%s' for write\n",filename);
   return(1);
   }

if(fwrite((void*)ds,1,1024,fp) != 1024) {
   fprintf(stderr,"ERROR: Incomplete write of header in write_datset()\n");
   fclose(fp);
   return(1);
   }

data_size = ds->nx * ds->ny * ds->nz;

switch(ds->mode) {
 case WU_BYTE: break;
 case WU_FLOAT:  data_size *= sizeof(float);   break;
 case WU_USHORT: data_size *= sizeof(unsigned short); break;
 case WU_SHORT:  data_size *= sizeof(short);   break;
 case WU_INT:    data_size *= sizeof(int);     break;
 default:
   fprintf(stderr,"ERROR: Invalid mode in write_dataset()\n");
   fclose(fp);
   return(1);
 }

if(fwrite((void*)ds->data,1,data_size,fp) != data_size) {
  fprintf(stderr,"WARNING: Writing data for file `%s' incomplete\n",filename);
  }

if((ds->user16_footer_size > 0)&&(ds->footer != (char*)NULL)) 
 if(fwrite((void*)ds->footer,1,
   ds->user16_footer_size,fp) != ds->user16_footer_size) {
   fprintf(stderr,
        "WARNING: Writing footer for file `%s' incomplete\n",filename);
  }
fclose(fp);
return(0);
}
