/***************************************************************************
  COPYRIGHT 1996 BIOMEDICAL COMPUTER LABORATORY, WASHINGTON UNIVERSITY, MO
****************************************************************************/
#define CheckNull(dat) if(dat == NULL) { \
	fprintf(stderr,"ERROR: Out of memory.\n"); exit(1); }

/*
 Data types for 'mode' field
*/
#define WU_BYTE 0
#define	WU_SHORT 1
#define WU_FLOAT 2
#define WU_USHORT 4
#define	WU_INT 5

#define WU_MAXDIM	10000

/*
   Magic number in header to make sure header VALID 
   and byte swapped correctly
*/
#define TVAL	1234

/* What TVAL looks like if byte-swapped */
#define SWABBED_TVAL -771489792

/* used to manipulate the dirty bit */

#define STAT_FLAG_SET	0x00000001	
#define STAT_FLAG_RESET	0xFFFFFFFE	

#define headtype	int

/*   
     The osm_ds is a structure which contains all the elements found in an
     image file.  The first section of the structure is 1024 bytes.  This 
     corresponds to the MRC HEADER used by UCSF algorithms.  The remaining 
     fields allow the image file data and footer portions to be included in the
     OSM dataset structure
*/

typedef struct _wuhdr
{

/******     MRC_HEADER FORMAT BEGINS HERE    *****/

/****** start fields used by WASHU programs ******/

   headtype 		nx;		/* number of columns */
   headtype 		ny;		/* number of rows */
   headtype 		nz;		/* number of sections */

   headtype 		mode;		/* data type */

/****** end   fields used by WASHU programs ******/

   headtype 		xstart;		/* number of first column in map */
   headtype 		ystart;		/* number of first row in map */
   headtype 		zstart;		/* number of first section in map */

   headtype 		mx;		/* number of intervals along X */
   headtype 		my;		/* number of intervals along Y */
   headtype 		mz;		/* number of intervals along Z */


   float		xlength;	/* cell dimensions (cell/mxyz = */
   float		ylength;	/* pixel spacing)               */
   float		zlength; 

   float		alpha;		/* cell angles (degrees) */
   float		beta;
   float		gamma;

   headtype 		col_axis;	/* which axis corresponds to  */
   headtype 		row_axis;	/* columns, rows, sections    */
   headtype 		sect_axis;      /* (X=1,Y=2,Z=3)              */

/****** start fields used by WASHU programs ******/

   float		amin;		/* minimum intensity value */
   float		amax;		/* maximum intensity value */
   float		amean;		/* mean intensity value */

/****** end   fields used by WASHU programs ******/

   headtype 		ispg;		/* space group number */
   headtype 		nsymbt;		/* number of bytes used for symmetry */
					/* operators                         */
 
/****** start fields used by WASHU programs ******/

   headtype 		user1_flags;	/* used for dataset dirty bit */ 
   headtype 		user2;
   headtype 		user3;
   headtype 		user4;
   headtype 		user5;
   headtype 		user6;
   headtype 		user7;
   headtype 		user8;
   headtype 		user9;
   headtype 		user10;
   headtype 		user11;
   headtype 		user12;
   headtype 		user13;
   headtype 		user14;
   headtype 		user15;
   headtype 		user16_footer_size; /* size of the footer in bytes*/

/****** end   fields used by WASHU programs ******/

   short 		id_type;	/* type of data set */
   short 		lens_type;	/* lens type */

   short		n1;		/* unknown */
   short		n2;		/* unknown */
   headtype 		data_value2;	/* not a readable format */

   float		tilt1;		/* two tilt sets (original 1-3 and */
   float		tilt2;		/* current 3-6)                    */
   float		tilt3;
   float		tilt4;
   float		tilt5;
   float		tilt6;

   headtype 		wave1;		/* multiple wave length information */
   headtype 		wave2;
   headtype 		wave3;

   float		xorigin;	/* x,y,z origin of the image */
   float		yorigin;
   float		zorigin;

/****** start fields used by WASHU programs ******/

   headtype 		num_labels;	/* number of text strings being used */
   char			label[10][80];	/* storage for first ten strings */

/****** MRC_HEADER FORMAT ENDS HERE (1024 bytes) *****/

   char			*data;		/* pointer to the data */
   char			*footer;	/* pointer to the footer */

/****** end   fields used by WASHU programs ******/
} osm_ds;

/* header procedure prototypes */

extern void print_notes (osm_ds* ds, int start, int end);
/* where strlen (text) <= 79 */
extern void add_history (osm_ds* ds, char* text);	
extern void clear_notes (osm_ds* ds);
extern void update_stats (osm_ds* ds);
extern osm_ds *read_dataset(char* filename);
extern osm_ds *read_header(char *fname);
extern int write_dataset (osm_ds* ds, char* filename);
extern void add_comment (osm_ds* ds, char* text);
extern void get_comment (osm_ds* ds);
extern int tobyte(osm_ds *ds);
extern int tofloat(osm_ds *ds);
extern int toushort(osm_ds *ds);
extern void SwabInt(int *d);
