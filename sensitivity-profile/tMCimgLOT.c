/********************************************************************************
*          Monte-Carlo Simulation for Light Transport in 3D Volumes             *
*********************************************************************************
*                                                                               *
* Copyright (C) 2002-2008,  David Boas    (dboas <at> nmr.mgh.harvard.edu)      *
*               2008        Jay Dubb      (jdubb <at> nmr.mgh.harvard.edu)      *
*               2008        Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)      *
*               2017        L. Nathan Perkins (lnp <at> bu.edu)                 *
*                                                                               *
* License:  4-clause BSD License, see LICENSE for details                       *
*                                                                               *
* Example:                                                                      *
*         tMCimgLOT input                                                       *
********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define INPUT_VERSION 1

#define C_VACUUM 2.9979e11f
#define TRUE 1
#define FALSE 0
#define MIN(a,b) ((a)<(b)?(a):(b))
#define FP_DIV_ERR  1e-8f

/* MACRO FOR RANDOM 0 to 1 FLOAT */
#define RANDF() ((REAL)rand() / RAND_MAX)

/* MACRO TO CONVERT 3D INDEX TO LINEAR INDEX. */
#define mult2linear(i,j,k,a1,a3)  (((k)-Izmin)*nIxy+((j)-Iymin)*nIx+((i)-Ixmin)+a3*nIxyz+a1*nIxyza3)

#ifdef VOXSIZE_EQ_1
#define DIST2VOX(x,s) (((int)(x)))
#else
#define DIST2VOX(x,s) ((int)((x)*(s)))
#endif

#define EPS 2.2204e-16f

#ifdef SINGLE_PREC
#define READ_ONE_REAL "%f"
#define READ_THREE_REALS "%f %f %f"
#define READ_REAL_INT_INT_INT "%f %d %d %d"
typedef float REAL;
#else
#define READ_ONE_REAL "%lf"
#define READ_THREE_REALS "%lf %lf %lf"
#define READ_REAL_INT_INT_INT "%lf %d %d %d"
typedef double REAL;
#endif

#define MAX_TISS_NUM  100
#define ASSERT(exp) tmc_assert(exp,__FILE__,__LINE__);

int idum; /* SEED FOR RANDOM NUMBER GENERATOR - A LARGE NEGATIVE NUMBER IS REQUIRED */
void tmc_error(int id, const char *msg, const char *fname, const int linenum);
void tmc_assert(int ret, const char *fname, const int linenum);

#ifdef ANGLE_FROM_NA
void tmc_perturb_angle(REAL *cx, REAL *cy, REAL *cz, REAL cna);
static int fiberNA(REAL NA, REAL sinth);
#endif

#define MAX_FILE_PATH 1024

int main(int argc, char *argv[])
{
    int version;
    
    int i, j, k, ii, jj, a1, a3;
    int N; /* NUMBER OF PHOTONS RUN SO FAR */
    int NT; /* TOTAL NUMBER OF PHOTONS TO RUN */
    int Ntissue; /* NUMBER OF TISSUE TYPES DESCRIBED IN THE IMAGE FILE */
    
    REAL foo, foo2; /* TEMPORARY VARIABLES */
    REAL ffoo;
    
    char ***tissueType; /* STORE THE IMAGE FILE */
    short tissueIndex;
    int nxstep, nystep, nzstep; /* DIMENSIONS OF THE IMAGE FILE */
    REAL xstep, ystep, zstep, rxstep, rystep, rzstep, minstepsize; /* VOXEL DIMENSIONS */
    int nA1step, nA3step;
    
    REAL tmus[MAX_TISS_NUM], tmua[MAX_TISS_NUM]; /* OPTICAL PROPERTIES OF THE DIFFERENT TISSUE TYPES */
    REAL tg[MAX_TISS_NUM], tn[MAX_TISS_NUM];
    
    REAL x,y,z; /* CURRENT PHOTON POSITION */
    REAL xi, yi, zi; /* INITIAL POSITION OF THE PHOTON */
    
    REAL gg, phi, theta, sphi, cphi, stheta, ctheta; /* SCATTERING ANGLES */
    REAL c1,c2,c3; /* DIRECTION COSINES */
    REAL c1o, c2o, c3o; /* OLD DIRECTION COSINES */
    REAL cxi, cyi, czi; /* INITIAL DIRECTION COSINES */
#ifdef ANGLE_FROM_NA
    REAL cna; /* NUMERICAL APERTURE FOR CALCULATING ANGLE */
#endif
    
    REAL *II, IIout[2]; /* FOR STORING THE 2-PT FLUENCE, IIout is for outside the II range */
    
    int Ixmin, Ixmax, Iymin, Iymax, Izmin, Izmax; /* MIN and MAX X,Y,Z FOR STORING THE 2-PT FLUENCE */
    int nIxstep, nIystep, nIzstep;
    int nIxyz,nIxy,nIx, nIxyza3, nIxyza13;
    
    REAL minT, maxT; /* MIN AND MAX TIME FOR SAMPLING THE 2-PT FLUENCE */
    REAL stepT, stepL; /* TIME STEP AND CORRESPONDING LENGTH STEP FOR SAMPLING THE 2-PT FLUENCE */
    REAL stepT_r, stepT_too_small; /* STEPT_R REMAINDER GATE WIDTH */
    
    REAL Lresid, Ltot, Lmin, Lmax, Lnext, step, nTstep_float;
    int nTstep, nTstep_int, tindex;
    
    int nDets; /* SPECIFY NUMBER OF DETECTORS*/
    REAL detRad; /* SPECIFY DETECTOR RADIUS */
    int **detLoc; /* and DETECTOR X,Y,Z LOCATIONS */
    REAL **detPos; /* and DETECTOR X,Y,Z LOCATIONS */
    
    
    REAL P2pt; /* PHOTON WEIGHT */
    
    REAL lenTiss[MAX_TISS_NUM]; /* THE LENGTH SPENT IN EACH TISSUE TYPE BY THE CURRENT PHOTON */
#ifdef MOMENTUM_TRANSFER
    REAL momTiss[MAX_TISS_NUM];
#endif
    REAL rnm; /* RANDOM NUMBER */
    
    FILE *fp; /* FILE POINTERS FOR SAVING THE DATA */
    char filenm[MAX_FILE_PATH]; /* FILE NAME FOR DATA FILE */
    char segFile[MAX_FILE_PATH]; /* FILE NAME FOR IMAGE FILE */
    
    int sizeof_lenTissArray;
#ifdef MOMENTUM_TRANSFER
    int sizeof_momTissArray;
#endif
    
    /* GET THE COMMAND LINE ARGUMENTS */
    if (argc != 2) {
        printf("usage: tMCimgLOT input_file (.inp assumed)\n");
        exit(1);
    }
    
    /*********************************************************
     OPEN AND READ THE INPUT FILE
     *********************************************************/
    sprintf(filenm, "%s.inp", argv[1]);
    printf("Loading configurations from %s\n", filenm);
    if ((fp = fopen(filenm, "r")) == NULL) {
        printf("usage: tMCimgLOT input_file (.inp assumed)\n");
        printf("input_file = %s does not exist.\n", filenm);
        exit(1);
    }
    
    /* READ THE INPUT FILE */
    ASSERT(fscanf(fp, "%d", &version) != 1); /* INPUT VERSION FORMAT */
    
    // check version
    if (version != INPUT_VERSION) {
        printf("Unexpected input file version (%d). Expected %d.\n", version, INPUT_VERSION);
        fclose(fp);
        exit(1);
    }
    
    ASSERT(fscanf(fp, "%d", &NT) != 1); /* TOTAL NUMBER OF PHOTONS */
    ASSERT(fscanf(fp, "%d", &idum) != 1); /* RANDOM NUMBER SEED */
    ASSERT(fscanf(fp, READ_THREE_REALS, &xi, &yi, &zi) != 3); /* INITIAL POSITION OF PHOTON */
    ASSERT(fscanf(fp, READ_THREE_REALS, &cxi, &cyi, &czi) != 3); /* INITIAL DIRECTION OF PHOTON */
#ifdef ANGLE_FROM_NA
    ASSERT(fscanf(fp, READ_ONE_REAL, &cna) != 1); /* NUMERICAL APERTURE OF PHOTON */
#endif
    ASSERT(fscanf(fp, READ_THREE_REALS, &minT, &maxT, &stepT) != 3); /* MIN, MAX, STEP TIME FOR RECORDING */
    ASSERT(fscanf(fp, "%d %d", &nA1step, &nA3step) != 2); /* NUMBER OF ANGULAR STEPS FOR II */
    
    /* Calculate number of gates, taking into account floating point division errors. */
    nTstep_float = (maxT - minT) / stepT;
    nTstep_int = (int)(nTstep_float);
    stepT_r = fabs(nTstep_float - nTstep_int) * stepT;
    stepT_too_small = FP_DIV_ERR * stepT;
    if (stepT_r < stepT_too_small)
        nTstep = nTstep_int;
    else
        nTstep = ceil(nTstep_float);
    
    ASSERT(fscanf(fp, "%s", segFile) != 1); /* FILE CONTAINING TISSUE STRUCTURE */
    
    /* READ IMAGE DIMENSIONS */
    ASSERT(fscanf(fp, READ_REAL_INT_INT_INT, &xstep, &nxstep, &Ixmin, &Ixmax) != 4);
    ASSERT(fscanf(fp, READ_REAL_INT_INT_INT, &ystep, &nystep, &Iymin, &Iymax) != 4);
    ASSERT(fscanf(fp, READ_REAL_INT_INT_INT, &zstep, &nzstep, &Izmin, &Izmax) != 4);
    --Ixmin; --Ixmax; --Iymin; --Iymax; --Izmin; --Izmax;
    nIxstep = Ixmax - Ixmin + 1;
    nIystep = Iymax - Iymin + 1;
    nIzstep = Izmax - Izmin + 1;
    
    minstepsize = MIN(xstep, MIN(ystep, zstep)); /* get the minimum dimension */
    rxstep = 1.f / xstep;
    rystep = 1.f / ystep;
    rzstep = 1.f / zstep;
    
    /* seed random number generator */
    if (idum != 0) {
        srand(abs(idum));
    } else {
        idum = time(NULL);
        srand(idum);
    }
    
    /* READ NUMBER OF TISSUE TYPES AND THEIR OPTICAL PROPERTIES */
    ASSERT(fscanf(fp, "%d", &Ntissue) != 1);
    tmus[0] = -999.f; tmua[0] = -999.f; tg[0] = -999.f; tn[0] = -999.f;
    for (i = 1; i <= Ntissue; ++i) {
#ifdef SINGLE_PREC
        ASSERT(fscanf(fp, "%f %f %f %f", &tmus[i], &tg[i], &tmua[i], &tn[i]) != 4);
#else
        ASSERT(fscanf(fp, "%lf %lf %lf %lf", &tmus[i], &tg[i], &tmua[i], &tn[i]) != 4);
#endif
        if (fabs(tn[i] - 1.0f) > EPS) {
            printf("WARNING: The code does not yet support n != 1.0\n");
        }
        if (fabs(tmus[i]) < EPS) {
            printf("ERROR: The code does support mus = 0.0\n");
            exit(1);
        }
    }
    
    /* READ NUMBER OF DETECTORS, DETECTOR RADIUS, AND DETECTOR LOCATIONS */
#ifdef SINGLE_PREC
    ASSERT(fscanf(fp, "%d %f", &nDets, &detRad) != 2);
#else
    ASSERT(fscanf(fp, "%d %lf", &nDets, &detRad) != 2);
#endif
    detLoc = (int **)malloc(nDets * sizeof(int *));
    detPos = (REAL **)malloc(nDets * sizeof(REAL *));
    
    for (i = 0; i < nDets; ++i){
        detLoc[i] = (int *)malloc(3 * sizeof(int));
        detPos[i] = (REAL *)malloc(3 * sizeof(REAL));
    }
    for (i = 0; i < nDets; ++i) {
        ASSERT(fscanf(fp, READ_THREE_REALS, detPos[i], detPos[i] + 1, detPos[i] + 2) != 3);
        detLoc[i][0] = (int)(detPos[i][0] * rxstep) - 1;
        detLoc[i][1] = (int)(detPos[i][1] * rystep) - 1;
        detLoc[i][2] = (int)(detPos[i][2] * rzstep) - 1;
    }
    
    fclose(fp);
    
    /* NORMALIZE THE DIRECTION COSINE OF THE SOURCE */
    foo = sqrt(cxi * cxi + cyi * cyi + czi * czi);  /*foo is the input */
    cxi /= foo;
    cyi /= foo;
    czi /= foo;
    
    /* CALCULATE THE MIN AND MAX PHOTON LENGTH FROM THE MIN AND MAX PROPAGATION TIMES */
    Lmax = maxT * C_VACUUM / tn[1];
    Lmin = minT * C_VACUUM / tn[1];
    stepL = stepT * C_VACUUM / tn[1];
    
    printf("Loading target medium volume from %s\n", segFile);
    
    /* READ IN THE SEGMENTED DATA FILE */
    fp = fopen(segFile, "rb");
    if (fp == NULL) {
        printf("ERROR: The binary image file %s was not found!\n", segFile);
        exit(1);
    }
    tissueType = (char ***)malloc(nxstep * sizeof(char **));
    for (i = 0; i < nxstep; ++i) {
        tissueType[i] = (char **)malloc(nystep * sizeof(char *));
        for (j = 0; j < nystep; ++j) {
            tissueType[i][j] = (char *)malloc(nzstep * sizeof(char));
        }
    }
    for (k = 0; k < nzstep; ++k) {
        for (j = 0; j < nystep; ++j) {
            for (i = 0; i < nxstep; ++i) {
                ASSERT(fscanf(fp, "%c", &tissueType[i][j][k]) != 1);
            }
        }
    }
    
    fclose(fp);
    
    
    nIxyz = nIzstep * nIxstep * nIystep;
    nIxy = nIxstep * nIystep;
    nIx = nIxstep;
    nIxyza3 = nIxyz * nA3step;
    nIxyza13 = nIxyz * nA1step * nA3step;
    
    /* ALLOCATE SPACE FOR AND INITIALIZE THE PHOTON FLUENCE TO 0 */
    II = (REAL *)malloc(nIxyza13 * nTstep * sizeof(REAL));
    memset((void*)II, 0, nIxyza13 * nTstep * sizeof(REAL));
    IIout[0] = 0.f;
    IIout[1] = 0.f;
    
    /* MAKE SURE THE SOURCE IS AT AN INTERFACE */
    i = DIST2VOX(xi, rxstep);
    j = DIST2VOX(yi, rystep);
    k = DIST2VOX(zi, rzstep);
    
    /* NUMBER PHOTONS EXECUTED SO FAR */
    N = 0;
    
    /* OPEN A FILE POINTER TO SAVE THE HISTORY INFORMATION */
    sprintf(filenm, "%s.his", argv[1]);
    fp = fopen(filenm, "wb");
    
    sizeof_lenTissArray = sizeof(REAL) * (Ntissue + 1);
#ifdef MOMENTUM_TRANSFER
    sizeof_momTissArray = sizeof(REAL) * (Ntissue + 1);
#endif
    
    /*********************************************************
     START MIGRATING THE PHOTONS
     GENERATING PHOTONS UNTIL NUMBER OF PHOTONS EXECUTED
     (N) IS EQUAL TO THE NUMBER TO BE GENERATED (NT)
     *********************************************************/
    
    printf("Launching %d photons\n", NT);
    
    while (N < NT) {
        ++N;
        
        /* SET THE PHOTON WEIGHT TO 1 AND INITIALIZE PHOTON LENGTH PARAMETERS */
        P2pt = 1.f;
        Ltot = 0.f;
        Lnext = minstepsize;
        Lresid = 0.f;
        
        /* INITIALIZE THE LENGTH IN EACH TISSUE TYPE */
        memset((void*)lenTiss, 0, sizeof_lenTissArray);
#ifdef MOMENTUM_TRANSFER
        memset((void*)momTiss, 0, sizeof_momTissArray);
#endif
        
        /* INITIAL SOURCE POSITION */
        x = xi;
        y = yi;
        z = zi;
        
        /* INITIAL DIRECTION OF PHOTON */
        c1 = cxi;
        c2 = cyi;
        c3 = czi;
        
#ifdef ANGLE_FROM_NA
        if (cna > 0) {
            /* LAUNCH WITHIN A SPECIFIC NA ALONG Z-AXIS */
            tmc_perturb_angle(&c1, &c2, &c3, cna);
            
            /* NORMALIZE THE DIRECTION COSINE OF THE SOURCE */
            // probably no longer needed, perturb seems to preserve the direction
            foo = sqrt(c1 * c1 + c2 * c2 + c3 * c3); /*foo is the input */
            c1 /= foo;
            c2 /= foo;
            c3 /= foo;
        }
#endif
        
        c1o = c1;
        c2o = c2;
        c3o = c3;
        
        /* PROPAGATE THE PHOTON */
        i = DIST2VOX(x, rxstep);
        j = DIST2VOX(y, rystep);
        k = DIST2VOX(z, rzstep);
        
        a3 = (int)round((float)nA3step * (c3 + 1.f) / 2.f - 0.5f);
        if (a3 == nA3step) {
            a3 = nA3step - 1;
        } else if (a3 < 0) {
            a3 = 0;
        }
        
        foo2 = atan2f(c2, c1);
        if (foo2 < 0.f) foo2 += 2.f * M_PI;
        a1 = (int)round((float)nA1step * foo2 / (2.f * M_PI) - 0.5f);
        if (a1 == nA1step) {
            a1 = nA1step - 1;
        } else if (a1 < 0) {
            printf("a1<0 : %d %.2f\n", a1, foo2);
            a1 = 0;
        }
        
        /* LOOP UNTIL TIME EXCEEDS GATE OR PHOTON ESCAPES */
        //        while ( Ltot<Lmax && i>=0 && i<nxstep && j>=0 && j<nystep && k>=0 && k<nzstep && (tissueIndex=tissueType[i][j][k])!=0 ) {
        while (Ltot < Lmax && k >= 1) {
            tissueIndex = 1;
            if (i >= 0 && i < nxstep && j >= 0 && j < nystep && k >= 0 && k < nzstep) {
                tissueIndex = tissueType[i][j][k];
            }
            
            /* CALCULATE SCATTERING LENGTH */
            rnm = RANDF(); /*ran( &idum, &ncall );*/
            if (rnm > EPS)
                Lresid = -log(rnm);
            else
                Lresid = -log(EPS);
            
            /* PROPAGATE THE PHOTON */
            //            while( Ltot<Lmax && Lresid>0. && i>=0 && i<nxstep && j>=0 && j<nystep && k>=0 && k<nzstep && (tissueIndex=tissueType[i][j][k])!=0 ) {
            while (Ltot < Lmax && Lresid > 0.f && k >= 1) {
                tissueIndex = 1;
                if (i >= 0 && i < nxstep && j >= 0 && j < nystep && k >= 0 && k < nzstep) {
                    tissueIndex = tissueType[i][j][k];
                }
                
                
                if (Ltot > Lnext && Ltot > Lmin) {
                    tindex = (int)((Ltot - Lmin) / stepL);
                    if (i >= Ixmin && i <= Ixmax && j >= Iymin && j <= Iymax && k >= Izmin && k <= Izmax && tindex < nTstep) {
#ifdef DEBUG
                        /*
                         printf("Scoring vox(%d,%d,%d) from photon %d at position (%0.1f,%0.1f,%0.1f) with fluence %f and direction (%0.1f,%0.1f,0.1%f)\n",
                         i, j, k, N, x, y, z, P2pt, c1, c2, c3);
                         */
#endif
                        II[mult2linear(i, j, k, a1, a3)] += P2pt;
                    } else {
                        IIout[0] += P2pt;
                    }
                    Lnext += minstepsize;
                }
                
                /*if scattering length is likely within a voxel, i.e. jump inside one voxel*/
                if ((foo = tmus[tissueIndex]) * minstepsize > Lresid) {
                    step = Lresid / foo;
                    x += c1 * step;
                    y += c2 * step;
                    z += c3 * step;
                    Ltot += step;
                    
                    P2pt *= exp(-tmua[tissueIndex] * step);
                    
                    lenTiss[tissueIndex] += (REAL)step;
                    Lresid = 0.f;
                } else {   /*if scattering length is bigger than a voxel, then move 1 voxel*/
                    x += c1 * minstepsize;
                    y += c2 * minstepsize;
                    z += c3 * minstepsize;
                    Ltot += minstepsize;
                    
                    P2pt *= exp(-tmua[tissueIndex] * minstepsize);
                    
                    Lresid -= foo * minstepsize;
                    lenTiss[tissueIndex] += minstepsize;
                }
                
                i = DIST2VOX(x, rxstep);
                j = DIST2VOX(y, rystep);
                k = DIST2VOX(z, rzstep);
                
            } /* PROPAGATE PHOTON */
            
            if (tissueIndex) {
                /* CALCULATE THE NEW SCATTERING ANGLE USING HENYEY-GREENSTEIN */
                gg = tg[tissueIndex];
                
                rnm = RANDF(); /*ran( &idum, &ncall );*/
                phi=2.0f * M_PI * rnm;
                cphi=cos(phi);
                sphi=sin(phi);
                
                rnm = RANDF(); /*ran( &idum, &ncall );*/
                if (gg > EPS) {
                    foo = (1.f - gg * gg) / (1.f - gg + 2.f * gg * rnm);
                    foo = foo * foo;
                    foo = (1.f + gg * gg - foo) / (2.f * gg);
                    theta = acos(foo);
                    stheta = sin(theta);
                    ctheta = foo;
                } else {  /*if g is exactly zero, then use isotropic scattering angle*/
                    theta = 2.0f * M_PI * rnm;
                    stheta = sin(theta);
                    ctheta = cos(theta);
                }
#ifdef MOMENTUM_TRANSFER
                if (theta > 0.f)
                    momTiss[tissueIndex] += 1.f - ctheta;
#endif
                c1o = c1;
                c2o = c2;
                c3o = c3;
                if (c3 < 1.f && c3 > -1.f) {
                    c1 = stheta * (c1o * c3o * cphi - c2o * sphi) / sqrt(1.f - c3o * c3o) + c1o * ctheta;
                    c2 = stheta * (c2o * c3o * cphi + c1o * sphi) / sqrt(1.f - c3o * c3o) + c2o * ctheta;
                    c3 = -stheta * cphi * sqrt(1 - c3o * c3o) + c3o * ctheta;
                } else {
                    c1 = stheta*cphi;
                    c2 = stheta*sphi;
                    c3 = ctheta*c3;
                }
                
                /* INDEX OF PHOTON DIRECTION */
                a3 = (int)round((float)nA3step * (c3 + 1.f) / 2.f - 0.5f );
                if (a3 == nA3step) {
                    a3 = nA3step - 1;
                } else if (a3 < 0) {
                    a3 = 0;
                }
                
                foo2 = atan2(c2, c1);
                if (foo2 < 0.f) foo2 += 2.f * M_PI;
                a1 = (int)round((float)nA1step * foo2 / (2.f * M_PI) - 0.5f);
                if (a1 == nA1step) {
                    a1 = nA1step - 1;
                } else if (a1 < 0) {
                    printf("a1<0 : %d %.2f\n", a1, foo2);
                    a1 = 0;
                }
                
            } /* LOOP UNTIL END OF SINGLE PHOTON */
        }
        
        /* SCORE EXITING PHOTON AND SAVE HISTORY FILES*/
        i = DIST2VOX(x, rxstep);
        j = DIST2VOX(y, rystep);
        k = DIST2VOX(z, rzstep);
        
        if (i >= 0 && i < nxstep && j >= 0 && j < nystep && k >= 0 && k < nzstep) {
            tissueIndex = tissueType[i][j][k];
            if (tissueIndex == 0) {
                tindex = (int)((Ltot - Lmin) / stepL);
                if (i >= Ixmin && i <= Ixmax && j >= Iymin && j <= Iymax && k >= Izmin && k <= Izmax && tindex < nTstep) {
#ifdef DEBUG
                    /*
                     printf("Scoring air vox(%d,%d,%d) from photon %d at position (%0.1f,%0.1f,%0.1f) with fluence %f and direction (%0.1f,%0.1f,0.1%f)\n",
                     i, j, k, N, x, y, z, P2pt, c1, c2, c3);
                     */
#endif
                    II[mult2linear(i,j,k,a1,a3)] -= P2pt;
                }
                
                /* LOOP THROUGH NUMBER OF DETECTORS */
                for (ii = 0; ii < nDets; ++ii) {
                    if (abs(i - detLoc[ii][0]) <= detRad) {
                        if (abs(j - detLoc[ii][1]) <= detRad) {
                            if (abs(k - detLoc[ii][2]) <= detRad) {
                                
                                /* WRITE TO THE HISTORY FILE */
                                ffoo = ii;
                                fwrite( &ffoo, sizeof(REAL), 1, fp );
                                for (jj = 1; jj <= Ntissue; ++jj) {
                                    fwrite(&lenTiss[jj], sizeof(REAL), 1, fp );
                                }
#ifdef MOMENTUM_TRANSFER
                                for (jj = 1; jj <= Ntissue; ++jj) {
                                    fwrite(&momTiss[jj], sizeof(REAL), 1, fp );
                                }
#endif
                            }
                        }
                    }
                }
                
                /* IF NO DETECTORS THEN SAVE EXIT POSITION */
                if (nDets == 0) {
                    ffoo = i; fwrite(&ffoo, sizeof(REAL), 1, fp);
                    ffoo = j; fwrite(&ffoo, sizeof(REAL), 1, fp);
                    ffoo = k; fwrite(&ffoo, sizeof(REAL), 1, fp);
                    for (jj = 1; jj <= Ntissue; ++jj) {
                        fwrite(&lenTiss[jj], sizeof(REAL), 1, fp);
                    }
#ifdef MOMENTUM_TRANSFER
                    for (jj = 1; jj <= Ntissue; ++jj) {
                        fwrite(&momTiss[jj], sizeof(REAL), 1, fp);
                    }
#endif
                }
                
            } /* End tissueIndex==0 */
        } else {
            IIout[1] -= P2pt;
        } /* End score exiting photon */
        
        
    } /* LOOP UNTIL ALL PHOTONS EXHAUSTED */
    
    /* CLOSE HISTORY FILE */
    fclose(fp); 
    
    
    /* SAVE FLUENCE DATA */
    sprintf(filenm, "%s.2pt", argv[1]);
    printf("Save photon fluence distribution to %s\n", filenm);
    
    fp = fopen(filenm, "wb");
    if (fp != NULL) {
        fwrite(II, sizeof(REAL), nIxyza13 * nTstep, fp);
        fwrite(IIout, sizeof(REAL), 2, fp );
        fclose(fp);
    } else {
        printf("ERROR: unable to save to %s\n", filenm);
        exit(1);
    }
    
    /*CLEAN UP*/
    for (i = 0; i < nDets; ++i) {
        free(detLoc[i]);
        free(detPos[i]);
    }
    free(detLoc);
    free(detPos);
    for (i = 0; i < nxstep; ++i) {
        for (j = 0; j < nystep; ++j) {
            free(tissueType[i][j]);
        }
        free(tissueType[i]);
    }
    free(tissueType);
    free(II);
    
    return 0;
}

#ifdef ANGLE_FROM_NA
/* 
 * BASED ON ORIGINAL: http://www.nmr.mgh.harvard.edu/DOT/resources/tmcimg/ but there
 * is probably a better way to calculate the following.
 */
void tmc_perturb_angle(REAL *cx, REAL *cy, REAL *cz, REAL cna) {
    REAL th0, dth;
    REAL ph0, dph;
    REAL p1, p2, p3;
    
    /* Original input direction characterized as a pair of rotation matrices characterized by two angles */
    
    th0 = acos(*cz);
    ph0 = atan2(*cy, *cx);
    
    /* Perturbation d\Omega to the propagation direction */
    
    dph = 2.f * M_PI * RANDF();
    
    do {
        dth = 2 * RANDF() - 1;
    } while (fiberNA(cna, dth) != 1);
    
    /* Convert to sin */
    dth = asin(dth);
    
    /* Convert perturbation into rectangular coordinates */
    p1 = sin(dth) * cos(dph);
    p2 = sin(dth) * sin(dph);
    p3 = cos(dth);
    
    *cx = cos(th0) * p1 + 0 * p2 + sin(th0) * p3;
    *cy = 0 * p1 + 1 * p2 + 0 * p3;
    *cz = -sin(th0) * p1 + 0 * p2 + cos(th0) * p3;
    
    p1 = *cx; /* Copy back to save a bit of storage space */
    p2 = *cy;
    p3 = *cz;
    
    *cx = cos(ph0) * p1 + sin(ph0) * p2 + 0 * p3;
    *cy = -sin(ph0) * p1 + cos(ph0) * p2 + 0 * p3;
    *cz = 0 * p1 + 0 * p2 + 1 * p3;
}

static int fiberNA(REAL NA, REAL sinth) {
    REAL a, b;
    
    if (sinth < -1 || sinth > 1) {
        ASSERT(1);
    }
    
    /* Delta-function distributions, will never match */
    if (NA <= 0.0 || NA > 1)
        return 0;
    
    /* Normalization for distribution */
    b = NA * NA / M_LN2;
    a = sqrt(b / M_PI) * erf(1 / sqrt(b));
    
    /* Accept if p < f(x) */
    
    if (RANDF() < exp(-b * sinth * sinth) / a)
        return 1;
    else
        return 0;
}
#endif

void tmc_error(int id,const char *msg,const char *fname,const int linenum) {
    fprintf(stderr, "tMCimg ERROR(%d):%s in %s:%d\n", id, msg, fname, linenum);
    exit(id);
}

void tmc_assert(int ret,const char *fname,const int linenum) {
    if(ret) tmc_error(ret, "assert error", fname, linenum);
}
