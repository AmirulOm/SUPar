/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SURANGE: $Revision: 1.19 $ ; $Date: 2013/06/24 16:56:54 $  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "thdr.h"
#include <omp.h>
#include <signal.h>

/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" SURANGE - get max and min values for non-zero header entries	",
" 								",
" surange <stdin	 					",
"								",
" Optional parameters:						",
"	key=		Header key(s) to range (default=all)	",
"	dim=0		dim seismic flag	",
"	    		0 = not dim, 1 = coord in ft, 2 = coord in m	",
" 								",
" Note: Gives partial results if interrupted			",
" 								",
" Output is: 							",
" number of traces 						",
" keyword min max (first - last) 				",
" north-south-east-west limits of shot, receiver and midpoint   ",
" if dim then also midpoint interval and line length   ",
    " 								",
NULL};

/* Credits:
 *  Arkansas: Chris Liner
 *              Added dim options Sept. 2017
 *  Stanford: Stewart A. Levin
 *              Added print of eastmost, northmost, westmost,
 *              southmost coordinates of shots, receivers, and 
 *              midpoints.  These coordinates have had any
 *              nonzero coscal header value applied.
 *	Geocon: Garry Perratt (output one header per line;
 *		option to specify headers to range;
 *		added first & last values where min<max)
 *	Based upon original by:
 *		SEP: Stew Levin
 *		CWP: Jack K. Cohen
 *
 * Note: the use of "signal" is inherited from BSD days and may
 *       break on some UNIXs.  It is dicey in that the responsibility
 *	 for program termination is lateraled back to the main.
 *
 */
/**************** end self doc ***********************************/


/* Prototypes */
//void printrange(traceHeader *tpmin, traceHeader *tpmax, traceHeader *tpfirst, traceHeader *tplast);
void printrange(segy *tpmin, segy *tpmax, segy *tpfirst, segy *tplast);
static void closeinput(void);
double** malloc2dFloat(int dim1, int dim2);
static segy tr, trmin, trmax, trfirst, trlast;

int
main(int argc, char **argv)
{
	//int ntr;			/* number of traces		*/
	int nkeys=0;			/* number of keywords to range	*/
	Value val;			/* value of current keyword	*/
	Value valmin;			/* smallest seen so far		*/
    	Value valmax;			/* largest seen so far		*/
    
    	Value vallast;
    	Value valfirst;
    	Value valmin_tmp;			/* smallest seen so far		*/
    	Value valmax_tmp;			/* largest seen so far		*/			

	cwp_String type;		/* data type of keyword		*/
	cwp_String key[SU_NKEYS];	/* array of keywords		*/
	int dim;			/* dim line with coords in ft (1) or m (2) */

        double eastShot[2], westShot[2], northShot[2], southShot[2];
        double eastRec[2], westRec[2], northRec[2], southRec[2];
        double eastCmp[2], westCmp[2], northCmp[2], southCmp[2];
        double dcoscal = 1.0;
        double sx, sy, gx, gy, mx, my;
        double mx1=0.0, my1=0.0;
        double mx2=0.0, my2=0.0, dm=0.0, dmin=0.0, dmax=0.0, davg=0.0;
        int coscal = 1;

	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);

	/* Get "key" value */
	if ((nkeys=countparval("key"))!=0) {
		getparstringarray("key",key);
	}
    
    /* get dim param... 0 = not dim */
    if (!getparint("dim", &dim)) dim = 0;
    
        checkpars();

	/* Zero out values of trmin and trmax */
	memset((void *) &trmin, 0, sizeof(segy));
	memset( (void *) &trmax, 0, sizeof(segy));
        northShot[0] = southShot[0] = eastShot[0] = westShot[0] = 0.0;
        northShot[1] = southShot[1] = eastShot[1] = westShot[1] = 0.0;
        northRec[0] = southRec[0] = eastRec[0] = westRec[0] = 0.0;
        northRec[1] = southRec[1] = eastRec[1] = westRec[1] = 0.0;
        northCmp[0] = southCmp[0] = eastCmp[0] = westCmp[0] = 0.0;
        northCmp[1] = southCmp[1] = eastCmp[1] = westCmp[1] = 0.0;
        sx = sy = gx = gy = mx = my = 0.0;

	/* Set up closing commands */
	signal(SIGINT, (void (*) (int)) closeinput);
    signal(SIGTERM, (void (*) (int)) closeinput);
    int totalSize = 0;
    int numSample = 0;
    
    if (!gettr(&tr)) err("can't get first trace");
    {
        numSample = tr.ns;
        int sz = fgetsizetr();
        totalSize = sz / ((numSample * 4) + 240);
    }
    
    thdr *trs = (thdr*) malloc(sizeof(thdr) * totalSize);
    memset(trs,0,sizeof(thdr)*totalSize);
    memcpy(&trs[0],&tr,sizeof(thdr));
    gettr(&tr);
    int ind=1;
    
    while (gettr(&tr)) {
	memcpy(&trs[ind],&tr,sizeof(thdr));
        ind++;
    }

    //Local variable for each thread.
    //TODO : set the number of thread less than num of trace
    const int numOfThread = omp_get_max_threads();
    omp_set_num_threads(numOfThread);
    printf("Numof thread = %i\n", numOfThread);

    segy *trmin_l = (segy*)malloc(sizeof(segy)*numOfThread);
    segy *trmax_l = (segy*)malloc(sizeof(segy)*numOfThread);
    segy *trfirst_l = (segy*)malloc(sizeof(segy)*numOfThread); 
    segy *trlast_l = (segy*)malloc(sizeof(segy)*numOfThread);
    memset(trmin_l,0,sizeof(segy)*numOfThread);
    memset(trmax_l,0,sizeof(segy)*numOfThread);
    memset(trfirst_l,0,sizeof(segy)*numOfThread);
    memset(trlast_l,0,sizeof(segy)*numOfThread);
    
    double **eastShot_l, **westShot_l, **northShot_l, **southShot_l;
    double **eastRec_l,  **westRec_l,  **northRec_l,  **southRec_l;
    double **eastCmp_l,  **westCmp_l,  **northCmp_l,  **southCmp_l;    
    double *dmin_l, *dmax_l, *davg_l;
    
    eastShot_l = malloc2dFloat(numOfThread,2); westShot_l = malloc2dFloat(numOfThread,2); northShot_l = malloc2dFloat(numOfThread,2); southShot_l = malloc2dFloat(numOfThread,2);
    eastRec_l = malloc2dFloat(numOfThread,2);  westRec_l = malloc2dFloat(numOfThread,2);  northRec_l = malloc2dFloat(numOfThread,2);  southRec_l = malloc2dFloat(numOfThread,2);
    eastCmp_l = malloc2dFloat(numOfThread,2);  westCmp_l = malloc2dFloat(numOfThread,2);  northCmp_l = malloc2dFloat(numOfThread,2);  southCmp_l = malloc2dFloat(numOfThread,2);
    dmin_l  = malloc(sizeof(double)*numOfThread); dmax_l = malloc(sizeof(double)*numOfThread); davg_l = malloc(sizeof(double)*numOfThread);

     #pragma omp parallel shared( trmin_l, trmax_l, trfirst_l, trlast_l, \
         eastShot_l, westShot_l, northShot_l, southShot_l, eastRec_l,  westRec_l,  \
         northRec_l,  southRec_l, eastCmp_l,  westCmp_l,  northCmp_l,  southCmp_l, dmin_l, dmax_l,trs) \
         private( dcoscal, coscal, sx, sy, gx, gy, mx, my, mx1, my1, mx2, my2, dm, davg,val,valmin,valmax,type) 
     {
         int tid = omp_get_thread_num();
         int lid = 0;
         dcoscal = 1.0;
         coscal = 1;
         mx1=0.0; my1=0.0;
         mx2=0.0; my2=0.0; dm=0.0; dmin=0.0, dmax=0.0, davg=0.0;
         sx = sy = gx = gy = mx = my = 0.0;
	     segy ltr;
         northShot_l[tid][0] = southShot_l[tid][0] = eastShot_l[tid][0] = westShot_l[tid][0] = 0.0;
         northShot_l[tid][1] = southShot_l[tid][1] = eastShot_l[tid][1] = westShot_l[tid][1] = 0.0;
         northRec_l[tid][0] = southRec_l[tid][  0] = eastRec_l[tid][0] = westRec_l[tid][0] = 0.0;
         northRec_l[tid][1] = southRec_l[tid][1] = eastRec_l[tid][1] = westRec_l[tid][1] = 0.0;
         northCmp_l[tid][0] = southCmp_l[tid][0] = eastCmp_l[tid][0] = westCmp_l[tid][0] = 0.0;
         northCmp_l[tid][1] = southCmp_l[tid][1] = eastCmp_l[tid][1] = westCmp_l[tid][1] = 0.0;

        //Each thread initialize first trace
        
          
         int ntr = 0;
         //For all of the traces,
         #pragma omp for schedule(static)
         for( lid=0; lid<totalSize ; lid ++)
         {
         register int i;
	  memcpy(&ltr, &trs[lid], sizeof(thdr)); 
         sx = sy = gx = gy = mx = my = 0.0;
	if(ntr==0){
         if (nkeys==0) {
             for (i = 0; i < SU_NKEYS; ++i) {
                 //gethval(&trs[tid], i, &val);
                 gethval(&ltr, i, &val);
                 puthval(&trmin_l[tid], i, &val);
                 puthval(&trmax_l[tid], i, &val);
                 puthval(&trfirst_l[tid], i, &val);
                                 if(i == 20) { coscal = val.h; if(coscal == 0) coscal = 1; dcoscal = (coscal > 0) ? 1.0*coscal : 1.0/coscal; }
                                 if(i == 21) sx = eastShot_l[tid][0] = westShot_l[tid][0] = northShot_l[tid][0] = southShot_l[tid][0] = val.i*dcoscal;
                                 if(i == 22) sy = eastShot_l[tid][1] = westShot_l[tid][1] = northShot_l[tid][1] = southShot_l[tid][1] = val.i*dcoscal;
                                 if(i == 23) gx = eastRec_l[tid][0] = westRec_l[tid][0] = northRec_l[tid][0] = southRec_l[tid][0] = val.i*dcoscal;
                                 if(i == 24) gy = eastRec_l[tid][1] = westRec_l[tid][1] = northRec_l[tid][1] = southRec_l[tid][1] = val.i*dcoscal;

             }
         } else	{
             register int j;
             for (i=0;i<nkeys;i++) {
                 j = getindex(key[i]);
                 //gethval(&trs[tid], j, &val);
                 gethval(&ltr, j, &val);
                 puthval(&trmin_l[tid], j, &val);
                 puthval(&trmax_l[tid], j, &val);
                 puthval(&trfirst_l[tid], j, &val);
             }
         }

         if(nkeys == 0) {
             mx = eastCmp_l[tid][0] = westCmp_l[tid][0] = northCmp_l[tid][0] = southCmp_l[tid][0] = 0.5*(eastShot_l[tid][0]+eastRec_l[tid][0]);
             my = eastCmp_l[tid][1] = westCmp_l[tid][1] = northCmp_l[tid][1] = southCmp_l[tid][1] = 0.5*(eastShot_l[tid][1]+eastRec_l[tid][1]);
         }
	}
	else{
           //  register int i;
           //          sx = sy = gx = gy = mx = my = 0.0;
	 //    memcpy(&ltr, &trs[lid], sizeof(thdr)); 
             if (nkeys==0) {
                     for (i = 0; i < SU_NKEYS; ++i) {
                     type = hdtype(getkey(i));
                     //gethval(&trs[lid], i, &val);
                     gethval(&ltr, i, &val);
                     gethval(&trmin_l[tid], i, &valmin);
                     gethval(&trmax_l[tid], i, &valmax);
                     if (valcmp(type, val, valmin) < 0)
			{ 
                         puthval(&trmin_l[tid], i, &val);
        		 if(tid == 0 && i == 0){
        		   printfval(type,val);
        		   printf("  %i\n",lid);
        		 }   
			}
                     if (valcmp(type, val, valmax) > 0)
                         puthval(&trmax_l[tid], i, &val);
                     puthval(&trlast_l[tid], i, &val);
                                     if(i == 20) { coscal = val.h; if(coscal == 0) coscal = 1; dcoscal = (coscal > 0) ? 1.0*coscal : 1.0/coscal; }
                                     if(i == 21)  sx = val.i*dcoscal;
                                     if(i == 22)  sy = val.i*dcoscal;
                                     if(i == 23)  gx = val.i*dcoscal;
                                     if(i == 24)  gy = val.i*dcoscal;
                 }
             } else	{
                 register int j;
                 for (i=0;i<nkeys;i++) {
                     type = hdtype(key[i]);
                     j = getindex(key[i]);
                     //gethval(&trs[lid], j, &val);
                     gethval(&ltr, j, &val);
                     gethval(&trmin_l[tid], j, &valmin);
                     gethval(&trmax_l[tid], j, &valmax);
                     if (valcmp(type, val, valmin) < 0)
                         puthval(&trmin_l[tid], j, &val);
                     if (valcmp(type, val, valmax) > 0)
                         puthval(&trmax_l[tid], j, &val);
                     puthval(&trlast_l[tid], j, &val);

                 }
             }

             if(nkeys == 0) {
                 mx = 0.5*(sx+gx); my = 0.5*(sy+gy);
                 if(eastShot_l[tid][0] < sx) {eastShot_l[tid][0] = sx; eastShot_l[tid][1] = sy;}
                 if(westShot_l[tid][0] > sx) {westShot_l[tid][0] = sx; westShot_l[tid][1] = sy;}
                 if(northShot_l[tid][1] < sy){northShot_l[tid][0] = sx; northShot_l[tid][1] = sy;}
                 if(southShot_l[tid][1] > sy){southShot_l[tid][0] = sx; southShot_l[tid][1] = sy;}
                 if(eastRec_l[tid][0] < gx) {eastRec_l[tid][0] = gx; eastRec_l[tid][1] = gy;}
                 if(westRec_l[tid][0] > gx) {westRec_l[tid][0] = gx; westRec_l[tid][1] = gy;}
                 if(northRec_l[tid][1] < gy){northRec_l[tid][0] = gx; northRec_l[tid][1] = gy;}
                 if(southRec_l[tid][1] > gy){southRec_l[tid][0] = gx; southRec_l[tid][1] = gy;}
                 if(eastCmp_l[tid][0] < mx) {eastCmp_l[tid][0] = mx; eastCmp_l[tid][1] = my;}
                 if(westCmp_l[tid][0] > mx) {westCmp_l[tid][0] = mx; westCmp_l[tid][1] = my;}
                 if(northCmp_l[tid][1] < my){northCmp_l[tid][0] = mx; northCmp_l[tid][1] = my;}
                 if(southCmp_l[tid][1] > my){southCmp_l[tid][0] = mx; southCmp_l[tid][1] = my;}
             }
          
             if (ntr == 1) {
                 /* get midpoint (mx1,my1) on trace 1 */
                 mx1 = 0.5*(ltr.sx+ltr.gx); 
                 my1 = 0.5*(ltr.sy+ltr.gy);
             }
             else if (ntr == 2) {
                 /* get midpoint (mx2,my2) on trace 2 */           
                 mx2 = 0.5*(ltr.sx+ltr.gx); 
                 my2 = 0.5*(ltr.sy+ltr.gy);
                 /* midpoint interval between traces 1 and 2 */
                 dm = sqrt( (mx1 - mx2)*(mx1 - mx2) + (my1 - my2)*(my1 - my2) );
                 /* set min, max and avg midpoint interval holders */
                 dmin_l[tid] = dm;
                 dmax_l[tid] = dm;
                 davg_l[tid] = (dmin_l[tid]+dmax_l[tid])/2.0;
                 /* hold this midpoint */
                 mx1 = mx2; 
                 my1 = my2;
             }
             else if (ntr > 2) {
                 /* get midpoint (mx,my) on this trace */           
                 mx2 = 0.5*(ltr.sx+ltr.gx); 
                 my2 = 0.5*(ltr.sy+ltr.gy);
                 /* get midpoint (mx,my) between this and previous trace */           
                 dm = sqrt( (mx1 - mx2)*(mx1 - mx2) + (my1 - my2)*(my1 - my2) );
                 /* reset min, max and avg midpoint interval holders, if needed */
                 if (dm < dmin_l[tid]) dmin_l[tid] = dm;
                 if (dm > dmax_l[tid]) dmax_l[tid] = dm;
                 davg_l[tid] = (davg_l[tid] + (dmin_l[tid]+dmax_l[tid])/2.0) / 2.0;
                 /* hold this midpoint */
                 mx1 = mx2; 
                 my1 = my2;
             }
          }

          ntr++;
         }
        
    }

    northShot[0] = northShot_l[0][0];
    southShot[0] = southShot_l[0][0];
    eastShot [0] = eastShot_l [0][0];
    westShot [0] = westShot_l [0][0];
    northShot[1] = northShot_l[0][1];
    southShot[1] = southShot_l[0][1];
    eastShot [1] = eastShot_l [0][1];
    westShot [1] = westShot_l [0][1];
    northRec [0] = northRec_l [0][0];
    southRec [0] = southRec_l [0][0];
    eastRec  [0] = eastRec_l  [0][0];
    westRec  [0] = westRec_l  [0][0];
    northRec [1] = northRec_l [0][1];
    southRec [1] = southRec_l [0][1];
    eastRec  [1] = eastRec_l  [0][1];
    westRec  [1] = westRec_l  [0][1];
    northCmp [0] = northCmp_l [0][0];
    southCmp [0] = southCmp_l [0][0];
    eastCmp  [0] = eastCmp_l  [0][0];
    westCmp  [0] = westCmp_l  [0][0];
    northCmp [1] = northCmp_l [0][1];
    southCmp [1] = southCmp_l [0][1];
    eastCmp  [1] = eastCmp_l  [0][1];
    westCmp  [1] = westCmp_l  [0][1];

    if (nkeys==0) {
        // can further optimize in just using memcpy
        for (int i = 0; i < SU_NKEYS; ++i) {
            type = hdtype(getkey(i));
            gethval(&trmin_l[0], i, &valmin);
            gethval(&trmax_l[0], i, &valmax);
            gethval(&trfirst_l[0], i, &valfirst);
            gethval(&trlast_l[0], i, &vallast);
            


            puthval(&trmin, i, &valmin);
            puthval(&trmax, i, &valmax);
            puthval(&trfirst, i, &valfirst);
            puthval(&trlast, i, &vallast); 

        }
    } else	{
        register int j;
        for (int i=0;i<nkeys;i++) {
            j = getindex(key[i]);

            gethval(&trmin_l[0], j, &valmin);
            gethval(&trmax_l[0], j, &valmax);
            gethval(&trfirst_l[0], j, &valfirst);
            gethval(&trlast_l[0], j, &vallast);

            puthval(&trmin, j, &valmin);
            puthval(&trmax, j, &valmax);
            puthval(&trfirst, j, &valfirst);
            puthval(&trlast, j, &vallast); 
        }
    }

    //Reduced all result to master thread
    for (int tid = 1; tid < numOfThread; tid++)
    {
        
        register int i;
        if (nkeys==0) {
                for (i = 0; i < SU_NKEYS; ++i) {
                type = hdtype(getkey(i));

                gethval(&trmin_l[tid], i, &valmin);
                gethval(&trmax_l[tid], i, &valmax);
                gethval(&trlast_l[tid], i, &vallast);
                
                //printfval(type,valmin);
                //printf("\n");

                gethval(&trmin, i, &valmin_tmp);
                gethval(&trmax, i, &valmax_tmp);

                if (valcmp(type, valmin, valmin_tmp) < 0)
                    puthval(&trmin, i, &valmin);
                if (valcmp(type, valmax, valmax_tmp) > 0)
                    puthval(&trmax, i, &valmax);
                puthval(&trlast, i, &vallast); 

            }
        } else	{
            register int j;
            for (i=0;i<nkeys;i++) {
                type = hdtype(key[i]);
                j = getindex(key[i]);

                gethval(&trmin_l[tid], j, &valmin);
                gethval(&trmax_l[tid], j, &valmax);
                gethval(&trlast_l[tid], j, &vallast);

                gethval(&trmin, j, &valmin_tmp);
                gethval(&trmax, j, &valmax_tmp);
                
                if (valcmp(type, valmin, valmin_tmp) < 0)
                    puthval(&trmin, j, &valmin);
                if (valcmp(type, valmax, valmax_tmp) > 0)
                    puthval(&trmax, j, &valmax);
                puthval(&trlast, j, &vallast); 

            }
        }

        if(nkeys == 0) {
            
            if(eastShot[0]  < eastShot_l[tid][0])  {eastShot[0]  = eastShot_l[tid][0];   eastShot[1]  = eastShot_l[tid][1];}
            if(westShot[0]  > westShot_l[tid][0])  {westShot[0]  = westShot_l[tid][0];   westShot[1]  = westShot_l[tid][1];}
            if(northShot[1] < northShot_l[tid][1]) {northShot[0] = northShot_l[tid][0];  northShot[1] = northShot_l[tid][1];}
            if(southShot[1] > southShot_l[tid][1]) {southShot[0] = southShot_l[tid][0];  southShot[1] = southShot_l[tid][1];}
            if(eastRec[0]   < eastRec_l[tid][0])   {eastRec[0]   = eastRec_l[tid][0];    eastRec[1]   = eastRec_l[tid][1];}
            if(westRec[0]   > westRec_l[tid][0])   {westRec[0]   = westRec_l[tid][0];    westRec[1]   = westRec_l[tid][1];}
            if(northRec[1]  < northRec_l[tid][1])  {northRec[0]  = northRec_l[tid][0];   northRec[1]  = northRec_l[tid][1] ;}
            if(southRec[1]  > southRec_l[tid][1])  {southRec[0]  = southRec_l[tid][0];   southRec[1]  = southRec_l[tid][1];}
            if(eastCmp[0]   < eastCmp_l[tid][0])   {eastCmp[0]   = eastCmp_l[tid][0];    eastCmp[1]   = eastCmp_l[tid][1];}
            if(westCmp[0]   > westCmp_l[tid][0])   {westCmp[0]   = westCmp_l[tid][0];    westCmp[1]   = westCmp_l[tid][1];}
            if(northCmp[1]  < northCmp_l[tid][1])  {northCmp[0]  = northCmp_l[tid][0];   northCmp[1]  = northCmp_l[tid][1];}
            if(southCmp[1]  > southCmp_l[tid][1])  {southCmp[0]  = southCmp_l[tid][0];   southCmp[1]  = southCmp_l[tid][1];}
        }

        if (dmin_l[tid] < dmin) dmin = dmin_l[tid];
        if (dmax_l[tid] > dmax) dmax = dmax_l[tid];
        //davg = (davg + (dmin_l[tid]+dmax_l[tid])/2.0) / 2.0;
        davg = (davg + davg_l[tid]) / 2.0; // <-- NOT SURE!!!

    }

    /* final davg is sum/elements
    davg = davg / (ntr-1); */
	int ntr = ind;
	printf("%d traces:\n",ntr);
	printrange(&trmin, &trmax, &trfirst, &trlast);
        if(nkeys == 0) {
            if(northShot[1] != 0.0 || southShot[1] != 0.0 ||
               eastShot[0] != 0.0 || westShot[0] != 0.0) printf(
                   "\nShot coordinate limits:\n"
                   "\tNorth(%g,%g) South(%g,%g) East(%g,%g) West(%g,%g)\n",
                   northShot[0],northShot[1],southShot[0],southShot[1],
                   eastShot[0],eastShot[1],westShot[0],westShot[1]);
            if(northRec[1] != 0.0 || southRec[1] != 0.0 ||
               eastRec[0] != 0.0 || westRec[0] != 0.0) printf(
                   "\nReceiver coordinate limits:\n"
                   "\tNorth(%g,%g) South(%g,%g) East(%g,%g) West(%g,%g)\n",
                   northRec[0],northRec[1],southRec[0],southRec[1],
                   eastRec[0],eastRec[1],westRec[0],westRec[1]);
            if(northCmp[1] != 0.0 || southCmp[1] != 0.0 ||
               eastCmp[0] != 0.0 || westCmp[0] != 0.0) printf(
                   "\nMidpoint coordinate limits:\n"
                   "\tNorth(%g,%g) South(%g,%g) East(%g,%g) West(%g,%g)\n",
                   northCmp[0],northCmp[1],southCmp[0],southCmp[1],
                   eastCmp[0],eastCmp[1],westCmp[0],westCmp[1]);
        }

    if (dim != 0){
        if (dim == 1) {
            printf("\n2D line: \n");
            printf("Min CMP interval = %g ft\n",dmin);
            printf("Max CMP interval = %g ft\n",dmax);
            printf("Line length = %g miles (using avg CMP interval of %g ft)\n",davg*ntr/5280,davg);
        }
        else if (dim == 2) {
            printf("ddim line: \n");
            printf("Min CMP interval = %g m\n",dmin);
            printf("Max CMP interval = %g m\n",dmax);
            printf("Line length = %g km (using avg CMP interval of %g m)\n",davg*ntr/1000,davg);
        }
    }
    
    return(CWP_Exit());
}

double** malloc2dFloat(int dim1, int dim2)
{
    double **t = (double**)malloc(sizeof(double*)*dim1);
    int i=0;
    for (i=0; i < dim1; i++)
    {
        t[i] = (double*)malloc(sizeof(double)*dim2);
    }

    return t;
}

/* printrange - print non-zero header values ranges	*/
void printrange(segy *tpmin, segy *tpmax, segy *tpfirst, segy *tplast)
{
	register int i = 0;
	Value valmin, valmax, valfirst, vallast;
	double dvalmin, dvalmax, dvalfirst, dvallast;
	cwp_String key;
	cwp_String type;
	int kmin = 0, kmax=SU_NKEYS;

	for (i = kmin; i < kmax; ++i) {
		key = getkey(i);
		type = hdtype(key);
		gethval(tpmin, i, &valmin);
		gethval(tpmax, i, &valmax);
		gethval(tpfirst, i, &valfirst);
		gethval(tplast, i, &vallast);
		dvalmin = vtod(type, valmin);
		dvalmax = vtod(type, valmax);
		dvalfirst = vtod(type, valfirst);
		dvallast = vtod(type, vallast);
		if (dvalmin || dvalmax) {
			if (dvalmin < dvalmax) {
				printf("%-8s ", key);
				printfval(type, valmin);
				printf(" ");
				printfval(type, valmax);
				printf(" (");
				printfval(type, valfirst);
				printf(" - ");
				printfval(type, vallast);
				printf(")");
			} else {
				printf("%-8s ", key);
				printfval(type, valmin);
			}
			putchar('\n');
		}
	}
	return;
}


static void closeinput(void) /* for graceful interrupt termination */
{
	/* Close stdin and open /dev/null in its place.  Now we are reading */
	/* from an empty file and the loops terminate in a normal fashion.  */

	efreopen("/dev/null", "r", stdin);
}
