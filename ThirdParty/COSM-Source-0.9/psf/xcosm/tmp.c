#define MAXPATHLEN 256

char lognm[MAXPATHLEN] = "psf_xcosm.log";
char psfnm[MAXPATHLEN] = "psf_xcosm.wu";
char prognm[MAXPATHLEN] = "CosmPsf";

float deltar = 0;
float deltaxy = 0;
float deltaxy_nyq = 0;

void WritePlane(int plane, float *dout, float *din, int w, int h, int indimx, int indimy)
{
int pls,i,j,deltx;
float *outloc,fval;
float *inloc;

deltx = indimx - w;
pls = w*h;
inloc = din;
outloc = dout + plane*pls;

for(j=0;j<h;j++,inloc += deltx)
 for(i=0;i<w;i++,inloc++,outloc++) {
        fval = *inloc;
        if(fval < 0.0) fval = 0.0;   /* force negative psf values to zero */
        *outloc = fval;
        }
}

void ZeroOut(float *f, int siz)
{
float *fptr;
int i;
for(fptr=f,i=0;i<siz;i++,fptr++)
  *fptr = 0.0;
}

int IsNotPowerOf2(int num)
{
int newnum,pcnt = 0;
newnum = num;
while(newnum > 0) { newnum >>= 1; pcnt++; }
pcnt--;
newnum = 1<<pcnt;
if(newnum != num) return(1);
else return(0);
}

