#include <unistd.h>
#include <stdio.h>
#include "idl_export.h"
    
#define NULL_VPTR ((IDL_VPTR) NULL)
    
#define GETVARADDR(v) ((v->flags&IDL_V_ARR) ? (void*)v->value.arr->data \
    : (void*) &v->value.c)
    
#define GETVARDATA(v,n) ((v->flags&IDL_V_ARR) \
    ? (*(n) = v->value.arr->n_elts, v->value.arr->data) \
    : (*(n) = 1, & v->value.c ) )

#define IDL_SINCE(maj,min) \
    ((IDL_VERSION_MAJOR>maj) || (IDL_VERSION_MAJOR==maj && IDL_VERSION_MINOR>=min))
    
/* 234567 */
/*   This fortran 77 file is to make interface between fortran 90 routines and IDL. */
extern  galaxy_idl_(void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*);

    
IDL_VPTR GALAXY_IDL(int argc, IDL_VPTR argv[])
{
  int in; /* Set flag appropriately before each variable block      */
  int i=0; /* Note it's use in indexing argv[i++]                   */
           /* This simplifies taking away extra arguments,          */
           /* typically those specifying the number of elements     */
           /* in input arrays (available as var->value.arr->n_elts) */

  int ndim;
  IDL_MEMINT dim[IDL_MAX_ARRAY_DIM];
  IDL_VPTR call[1]; /* For use when calling conversion routines */
  
  IDL_VPTR tmp;
/* ----- Parameter declarations ---------------------------------------- */
 /*  input/output variables                                               */

  IDL_VPTR NO_PHOTONS_=argv[i++];  /* NO_PHOTONS_ : REAL*8 :  :in? */
  IDL_VPTR NPRINT_=argv[i++];  /* NPRINT_ : INTEGER :  :in? */
  IDL_VPTR HGG_=argv[i++];  /* HGG_ : REAL*8 :  :in? */
  IDL_VPTR ALBEDO_=argv[i++];  /* ALBEDO_ : REAL*8 :  :in? */
  IDL_VPTR LUMINOSITY_=argv[i++];  /* LUMINOSITY_ : REAL*8 :  :in? */
  IDL_VPTR DUST1_=argv[i++];  /* DUST1_ : CHARACTER*(128) :  :in? */
  IDL_VPTR TAUFACE1_=argv[i++];  /* TAUFACE1_ : REAL*8 :  :in? */
  IDL_VPTR DUST_RSCALE1_=argv[i++];  /* DUST_RSCALE1_ : REAL*8 :  :in? */
  IDL_VPTR DUST_ZSCALE1_=argv[i++];  /* DUST_ZSCALE1_ : REAL*8 :  :in? */
  IDL_VPTR DUST_RMAX1_=argv[i++];  /* DUST_RMAX1_ : REAL*8 :  :in? */
  IDL_VPTR DUST_ZMAX1_=argv[i++];  /* DUST_ZMAX1_ : REAL*8 :  :in? */
  IDL_VPTR DUST2_=argv[i++];  /* DUST2_ : CHARACTER*(128) :  :in? */
  IDL_VPTR TAUFACE2_=argv[i++];  /* TAUFACE2_ : REAL*8 :  :in? */
  IDL_VPTR DUST_RSCALE2_=argv[i++];  /* DUST_RSCALE2_ : REAL*8 :  :in? */
  IDL_VPTR DUST_ZSCALE2_=argv[i++];  /* DUST_ZSCALE2_ : REAL*8 :  :in? */
  IDL_VPTR DUST_RMAX2_=argv[i++];  /* DUST_RMAX2_ : REAL*8 :  :in? */
  IDL_VPTR DUST_ZMAX2_=argv[i++];  /* DUST_ZMAX2_ : REAL*8 :  :in? */
  IDL_VPTR RMAX_=argv[i++];  /* RMAX_ : REAL*8 :  :in? */
  IDL_VPTR PMAX_=argv[i++];  /* PMAX_ : REAL*8 :  :in? */
  IDL_VPTR ZMAX_=argv[i++];  /* ZMAX_ : REAL*8 :  :in? */
  IDL_VPTR NR_=argv[i++];  /* NR_ : INTEGER :  :in? */
  IDL_VPTR NP_=argv[i++];  /* NP_ : INTEGER :  :in? */
  IDL_VPTR NZ_=argv[i++];  /* NZ_ : INTEGER :  :in? */
  IDL_VPTR DISK_NAME_=argv[i++];  /* DISK_NAME_ : CHARACTER*(128) :  :in? */
  IDL_VPTR DISK_RSCALE_=argv[i++];  /* DISK_RSCALE_ : REAL*8 :  :in? */
  IDL_VPTR DISK_ZSCALE_=argv[i++];  /* DISK_ZSCALE_ : REAL*8 :  :in? */
  IDL_VPTR DISK_RMAX_=argv[i++];  /* DISK_RMAX_ : REAL*8 :  :in? */
  IDL_VPTR DISK_ZMAX_=argv[i++];  /* DISK_ZMAX_ : REAL*8 :  :in? */
  IDL_VPTR BULGE_NAME_=argv[i++];  /* BULGE_NAME_ : CHARACTER*(128) :  :in? */
  IDL_VPTR SERSIC_INDEX_=argv[i++];  /* SERSIC_INDEX_ : REAL*8 :  :in? */
  IDL_VPTR REFF_=argv[i++];  /* REFF_ : REAL*8 :  :in? */
  IDL_VPTR AXIAL_RATIO_=argv[i++];  /* AXIAL_RATIO_ : REAL*8 :  :in? */
  IDL_VPTR BULGETODISK_=argv[i++];  /* BULGETODISK_ : REAL*8 :  :in? */
  IDL_VPTR INCLINATION_ANGLE_=argv[i++];  /* INCLINATION_ANGLE_ : REAL*8 :  :in? */
  IDL_VPTR POSITION_ANGLE_=argv[i++];  /* POSITION_ANGLE_ : REAL*8 :  :in? */
  IDL_VPTR PHASE_ANGLE_=argv[i++];  /* PHASE_ANGLE_ : REAL*8 :  :in? */
  IDL_VPTR DISTANCE_=argv[i++];  /* DISTANCE_ : REAL*8 :  :in? */
  IDL_VPTR NXIM_=argv[i++];  /* NXIM_ : INTEGER :  :in? */
  IDL_VPTR NYIM_=argv[i++];  /* NYIM_ : INTEGER :  :in? */
  IDL_VPTR DXIM_=argv[i++];  /* DXIM_ : REAL*8 :  :in? */
  IDL_VPTR DYIM_=argv[i++];  /* DYIM_ : REAL*8 :  :in? */
  IDL_VPTR IM_SCATT_=argv[i++];  /* IM_SCATT_ : REAL : (NXIM,NYIM) :in? */
  IDL_VPTR IM_DIREC_=argv[i++];  /* IM_DIREC_ : REAL : (NXIM,NYIM) :in? */
  IDL_VPTR IM_SCATT_SIG_=argv[i++];  /* IM_SCATT_SIG_ : REAL : (NXIM,NYIM) :in? */
  IDL_VPTR IM_DIREC_SIG_=argv[i++];  /* IM_DIREC_SIG_ : REAL : (NXIM,NYIM) :in? */
  IDL_VPTR PSF_FILE_=argv[i++];  /* PSF_FILE_ : CHARACTER*(128) :  :in? */
  IDL_VPTR LEFT_RIGHT_FOLD_=argv[i++];  /* LEFT_RIGHT_FOLD_ : INTEGER :  :in? */



  /* TYPE CHECKING / ALLOCATION SECTION */

  in = 1;           /* NO_PHOTONS_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(NO_PHOTONS_);
    IDL_ENSURE_SIMPLE(NO_PHOTONS_);
    IDL_EXCLUDE_STRING(NO_PHOTONS_);
    IDL_EXCLUDE_COMPLEX(NO_PHOTONS_);
    IDL_ENSURE_SCALAR(NO_PHOTONS_);
    call[0] = NO_PHOTONS_;
    NO_PHOTONS_ = IDL_CvtDbl(1,call); /* May cause NO_PHOTONS_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(NO_PHOTONS_); /* Output cannot be expression */
    IDL_StoreScalarZero(NO_PHOTONS_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* NPRINT_ : INTEGER :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(NPRINT_);
    IDL_ENSURE_SIMPLE(NPRINT_);
    IDL_EXCLUDE_STRING(NPRINT_);
    IDL_EXCLUDE_COMPLEX(NPRINT_);
    IDL_ENSURE_SCALAR(NPRINT_);
    call[0] = NPRINT_;
    NPRINT_ = IDL_CvtLng(1,call); /* May cause NPRINT_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(NPRINT_); /* Output cannot be expression */
    IDL_StoreScalarZero(NPRINT_,IDL_TYP_LONG);  /* Not for in/out! */
  }

  in = 1;           /* HGG_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(HGG_);
    IDL_ENSURE_SIMPLE(HGG_);
    IDL_EXCLUDE_STRING(HGG_);
    IDL_EXCLUDE_COMPLEX(HGG_);
    IDL_ENSURE_SCALAR(HGG_);
    call[0] = HGG_;
    HGG_ = IDL_CvtDbl(1,call); /* May cause HGG_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(HGG_); /* Output cannot be expression */
    IDL_StoreScalarZero(HGG_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* ALBEDO_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(ALBEDO_);
    IDL_ENSURE_SIMPLE(ALBEDO_);
    IDL_EXCLUDE_STRING(ALBEDO_);
    IDL_EXCLUDE_COMPLEX(ALBEDO_);
    IDL_ENSURE_SCALAR(ALBEDO_);
    call[0] = ALBEDO_;
    ALBEDO_ = IDL_CvtDbl(1,call); /* May cause ALBEDO_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(ALBEDO_); /* Output cannot be expression */
    IDL_StoreScalarZero(ALBEDO_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* LUMINOSITY_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(LUMINOSITY_);
    IDL_ENSURE_SIMPLE(LUMINOSITY_);
    IDL_EXCLUDE_STRING(LUMINOSITY_);
    IDL_EXCLUDE_COMPLEX(LUMINOSITY_);
    IDL_ENSURE_SCALAR(LUMINOSITY_);
    call[0] = LUMINOSITY_;
    LUMINOSITY_ = IDL_CvtDbl(1,call); /* May cause LUMINOSITY_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(LUMINOSITY_); /* Output cannot be expression */
    IDL_StoreScalarZero(LUMINOSITY_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DUST1_ : CHARACTER*(128) :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DUST1_);
    IDL_ENSURE_SIMPLE(DUST1_);
    IDL_EXCLUDE_COMPLEX(DUST1_);
    IDL_ENSURE_SCALAR(DUST1_);
    call[0] = DUST1_;
    DUST1_ = IDL_CvtByte(1,call); /* May cause DUST1_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DUST1_); /* Output cannot be expression */
    IDL_StoreScalarZero(DUST1_,IDL_TYP_BYTE);  /* Not for in/out! */
  }

  in = 1;           /* TAUFACE1_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(TAUFACE1_);
    IDL_ENSURE_SIMPLE(TAUFACE1_);
    IDL_EXCLUDE_STRING(TAUFACE1_);
    IDL_EXCLUDE_COMPLEX(TAUFACE1_);
    IDL_ENSURE_SCALAR(TAUFACE1_);
    call[0] = TAUFACE1_;
    TAUFACE1_ = IDL_CvtDbl(1,call); /* May cause TAUFACE1_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(TAUFACE1_); /* Output cannot be expression */
    IDL_StoreScalarZero(TAUFACE1_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DUST_RSCALE1_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DUST_RSCALE1_);
    IDL_ENSURE_SIMPLE(DUST_RSCALE1_);
    IDL_EXCLUDE_STRING(DUST_RSCALE1_);
    IDL_EXCLUDE_COMPLEX(DUST_RSCALE1_);
    IDL_ENSURE_SCALAR(DUST_RSCALE1_);
    call[0] = DUST_RSCALE1_;
    DUST_RSCALE1_ = IDL_CvtDbl(1,call); /* May cause DUST_RSCALE1_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DUST_RSCALE1_); /* Output cannot be expression */
    IDL_StoreScalarZero(DUST_RSCALE1_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DUST_ZSCALE1_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DUST_ZSCALE1_);
    IDL_ENSURE_SIMPLE(DUST_ZSCALE1_);
    IDL_EXCLUDE_STRING(DUST_ZSCALE1_);
    IDL_EXCLUDE_COMPLEX(DUST_ZSCALE1_);
    IDL_ENSURE_SCALAR(DUST_ZSCALE1_);
    call[0] = DUST_ZSCALE1_;
    DUST_ZSCALE1_ = IDL_CvtDbl(1,call); /* May cause DUST_ZSCALE1_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DUST_ZSCALE1_); /* Output cannot be expression */
    IDL_StoreScalarZero(DUST_ZSCALE1_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DUST_RMAX1_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DUST_RMAX1_);
    IDL_ENSURE_SIMPLE(DUST_RMAX1_);
    IDL_EXCLUDE_STRING(DUST_RMAX1_);
    IDL_EXCLUDE_COMPLEX(DUST_RMAX1_);
    IDL_ENSURE_SCALAR(DUST_RMAX1_);
    call[0] = DUST_RMAX1_;
    DUST_RMAX1_ = IDL_CvtDbl(1,call); /* May cause DUST_RMAX1_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DUST_RMAX1_); /* Output cannot be expression */
    IDL_StoreScalarZero(DUST_RMAX1_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DUST_ZMAX1_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DUST_ZMAX1_);
    IDL_ENSURE_SIMPLE(DUST_ZMAX1_);
    IDL_EXCLUDE_STRING(DUST_ZMAX1_);
    IDL_EXCLUDE_COMPLEX(DUST_ZMAX1_);
    IDL_ENSURE_SCALAR(DUST_ZMAX1_);
    call[0] = DUST_ZMAX1_;
    DUST_ZMAX1_ = IDL_CvtDbl(1,call); /* May cause DUST_ZMAX1_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DUST_ZMAX1_); /* Output cannot be expression */
    IDL_StoreScalarZero(DUST_ZMAX1_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DUST2_ : CHARACTER*(128) :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DUST2_);
    IDL_ENSURE_SIMPLE(DUST2_);
    IDL_EXCLUDE_COMPLEX(DUST2_);
    IDL_ENSURE_SCALAR(DUST2_);
    call[0] = DUST2_;
    DUST2_ = IDL_CvtByte(1,call); /* May cause DUST2_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DUST2_); /* Output cannot be expression */
    IDL_StoreScalarZero(DUST2_,IDL_TYP_BYTE);  /* Not for in/out! */
  }

  in = 1;           /* TAUFACE2_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(TAUFACE2_);
    IDL_ENSURE_SIMPLE(TAUFACE2_);
    IDL_EXCLUDE_STRING(TAUFACE2_);
    IDL_EXCLUDE_COMPLEX(TAUFACE2_);
    IDL_ENSURE_SCALAR(TAUFACE2_);
    call[0] = TAUFACE2_;
    TAUFACE2_ = IDL_CvtDbl(1,call); /* May cause TAUFACE2_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(TAUFACE2_); /* Output cannot be expression */
    IDL_StoreScalarZero(TAUFACE2_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DUST_RSCALE2_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DUST_RSCALE2_);
    IDL_ENSURE_SIMPLE(DUST_RSCALE2_);
    IDL_EXCLUDE_STRING(DUST_RSCALE2_);
    IDL_EXCLUDE_COMPLEX(DUST_RSCALE2_);
    IDL_ENSURE_SCALAR(DUST_RSCALE2_);
    call[0] = DUST_RSCALE2_;
    DUST_RSCALE2_ = IDL_CvtDbl(1,call); /* May cause DUST_RSCALE2_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DUST_RSCALE2_); /* Output cannot be expression */
    IDL_StoreScalarZero(DUST_RSCALE2_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DUST_ZSCALE2_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DUST_ZSCALE2_);
    IDL_ENSURE_SIMPLE(DUST_ZSCALE2_);
    IDL_EXCLUDE_STRING(DUST_ZSCALE2_);
    IDL_EXCLUDE_COMPLEX(DUST_ZSCALE2_);
    IDL_ENSURE_SCALAR(DUST_ZSCALE2_);
    call[0] = DUST_ZSCALE2_;
    DUST_ZSCALE2_ = IDL_CvtDbl(1,call); /* May cause DUST_ZSCALE2_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DUST_ZSCALE2_); /* Output cannot be expression */
    IDL_StoreScalarZero(DUST_ZSCALE2_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DUST_RMAX2_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DUST_RMAX2_);
    IDL_ENSURE_SIMPLE(DUST_RMAX2_);
    IDL_EXCLUDE_STRING(DUST_RMAX2_);
    IDL_EXCLUDE_COMPLEX(DUST_RMAX2_);
    IDL_ENSURE_SCALAR(DUST_RMAX2_);
    call[0] = DUST_RMAX2_;
    DUST_RMAX2_ = IDL_CvtDbl(1,call); /* May cause DUST_RMAX2_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DUST_RMAX2_); /* Output cannot be expression */
    IDL_StoreScalarZero(DUST_RMAX2_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DUST_ZMAX2_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DUST_ZMAX2_);
    IDL_ENSURE_SIMPLE(DUST_ZMAX2_);
    IDL_EXCLUDE_STRING(DUST_ZMAX2_);
    IDL_EXCLUDE_COMPLEX(DUST_ZMAX2_);
    IDL_ENSURE_SCALAR(DUST_ZMAX2_);
    call[0] = DUST_ZMAX2_;
    DUST_ZMAX2_ = IDL_CvtDbl(1,call); /* May cause DUST_ZMAX2_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DUST_ZMAX2_); /* Output cannot be expression */
    IDL_StoreScalarZero(DUST_ZMAX2_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* RMAX_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(RMAX_);
    IDL_ENSURE_SIMPLE(RMAX_);
    IDL_EXCLUDE_STRING(RMAX_);
    IDL_EXCLUDE_COMPLEX(RMAX_);
    IDL_ENSURE_SCALAR(RMAX_);
    call[0] = RMAX_;
    RMAX_ = IDL_CvtDbl(1,call); /* May cause RMAX_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(RMAX_); /* Output cannot be expression */
    IDL_StoreScalarZero(RMAX_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* PMAX_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(PMAX_);
    IDL_ENSURE_SIMPLE(PMAX_);
    IDL_EXCLUDE_STRING(PMAX_);
    IDL_EXCLUDE_COMPLEX(PMAX_);
    IDL_ENSURE_SCALAR(PMAX_);
    call[0] = PMAX_;
    PMAX_ = IDL_CvtDbl(1,call); /* May cause PMAX_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(PMAX_); /* Output cannot be expression */
    IDL_StoreScalarZero(PMAX_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* ZMAX_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(ZMAX_);
    IDL_ENSURE_SIMPLE(ZMAX_);
    IDL_EXCLUDE_STRING(ZMAX_);
    IDL_EXCLUDE_COMPLEX(ZMAX_);
    IDL_ENSURE_SCALAR(ZMAX_);
    call[0] = ZMAX_;
    ZMAX_ = IDL_CvtDbl(1,call); /* May cause ZMAX_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(ZMAX_); /* Output cannot be expression */
    IDL_StoreScalarZero(ZMAX_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* NR_ : INTEGER :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(NR_);
    IDL_ENSURE_SIMPLE(NR_);
    IDL_EXCLUDE_STRING(NR_);
    IDL_EXCLUDE_COMPLEX(NR_);
    IDL_ENSURE_SCALAR(NR_);
    call[0] = NR_;
    NR_ = IDL_CvtLng(1,call); /* May cause NR_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(NR_); /* Output cannot be expression */
    IDL_StoreScalarZero(NR_,IDL_TYP_LONG);  /* Not for in/out! */
  }

  in = 1;           /* NP_ : INTEGER :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(NP_);
    IDL_ENSURE_SIMPLE(NP_);
    IDL_EXCLUDE_STRING(NP_);
    IDL_EXCLUDE_COMPLEX(NP_);
    IDL_ENSURE_SCALAR(NP_);
    call[0] = NP_;
    NP_ = IDL_CvtLng(1,call); /* May cause NP_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(NP_); /* Output cannot be expression */
    IDL_StoreScalarZero(NP_,IDL_TYP_LONG);  /* Not for in/out! */
  }

  in = 1;           /* NZ_ : INTEGER :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(NZ_);
    IDL_ENSURE_SIMPLE(NZ_);
    IDL_EXCLUDE_STRING(NZ_);
    IDL_EXCLUDE_COMPLEX(NZ_);
    IDL_ENSURE_SCALAR(NZ_);
    call[0] = NZ_;
    NZ_ = IDL_CvtLng(1,call); /* May cause NZ_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(NZ_); /* Output cannot be expression */
    IDL_StoreScalarZero(NZ_,IDL_TYP_LONG);  /* Not for in/out! */
  }

  in = 1;           /* DISK_NAME_ : CHARACTER*(128) :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DISK_NAME_);
    IDL_ENSURE_SIMPLE(DISK_NAME_);
    IDL_EXCLUDE_COMPLEX(DISK_NAME_);
    IDL_ENSURE_SCALAR(DISK_NAME_);
    call[0] = DISK_NAME_;
    DISK_NAME_ = IDL_CvtByte(1,call); /* May cause DISK_NAME_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DISK_NAME_); /* Output cannot be expression */
    IDL_StoreScalarZero(DISK_NAME_,IDL_TYP_BYTE);  /* Not for in/out! */
  }

  in = 1;           /* DISK_RSCALE_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DISK_RSCALE_);
    IDL_ENSURE_SIMPLE(DISK_RSCALE_);
    IDL_EXCLUDE_STRING(DISK_RSCALE_);
    IDL_EXCLUDE_COMPLEX(DISK_RSCALE_);
    IDL_ENSURE_SCALAR(DISK_RSCALE_);
    call[0] = DISK_RSCALE_;
    DISK_RSCALE_ = IDL_CvtDbl(1,call); /* May cause DISK_RSCALE_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DISK_RSCALE_); /* Output cannot be expression */
    IDL_StoreScalarZero(DISK_RSCALE_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DISK_ZSCALE_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DISK_ZSCALE_);
    IDL_ENSURE_SIMPLE(DISK_ZSCALE_);
    IDL_EXCLUDE_STRING(DISK_ZSCALE_);
    IDL_EXCLUDE_COMPLEX(DISK_ZSCALE_);
    IDL_ENSURE_SCALAR(DISK_ZSCALE_);
    call[0] = DISK_ZSCALE_;
    DISK_ZSCALE_ = IDL_CvtDbl(1,call); /* May cause DISK_ZSCALE_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DISK_ZSCALE_); /* Output cannot be expression */
    IDL_StoreScalarZero(DISK_ZSCALE_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DISK_RMAX_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DISK_RMAX_);
    IDL_ENSURE_SIMPLE(DISK_RMAX_);
    IDL_EXCLUDE_STRING(DISK_RMAX_);
    IDL_EXCLUDE_COMPLEX(DISK_RMAX_);
    IDL_ENSURE_SCALAR(DISK_RMAX_);
    call[0] = DISK_RMAX_;
    DISK_RMAX_ = IDL_CvtDbl(1,call); /* May cause DISK_RMAX_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DISK_RMAX_); /* Output cannot be expression */
    IDL_StoreScalarZero(DISK_RMAX_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DISK_ZMAX_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DISK_ZMAX_);
    IDL_ENSURE_SIMPLE(DISK_ZMAX_);
    IDL_EXCLUDE_STRING(DISK_ZMAX_);
    IDL_EXCLUDE_COMPLEX(DISK_ZMAX_);
    IDL_ENSURE_SCALAR(DISK_ZMAX_);
    call[0] = DISK_ZMAX_;
    DISK_ZMAX_ = IDL_CvtDbl(1,call); /* May cause DISK_ZMAX_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DISK_ZMAX_); /* Output cannot be expression */
    IDL_StoreScalarZero(DISK_ZMAX_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* BULGE_NAME_ : CHARACTER*(128) :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(BULGE_NAME_);
    IDL_ENSURE_SIMPLE(BULGE_NAME_);
    IDL_EXCLUDE_COMPLEX(BULGE_NAME_);
    IDL_ENSURE_SCALAR(BULGE_NAME_);
    call[0] = BULGE_NAME_;
    BULGE_NAME_ = IDL_CvtByte(1,call); /* May cause BULGE_NAME_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(BULGE_NAME_); /* Output cannot be expression */
    IDL_StoreScalarZero(BULGE_NAME_,IDL_TYP_BYTE);  /* Not for in/out! */
  }

  in = 1;           /* SERSIC_INDEX_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(SERSIC_INDEX_);
    IDL_ENSURE_SIMPLE(SERSIC_INDEX_);
    IDL_EXCLUDE_STRING(SERSIC_INDEX_);
    IDL_EXCLUDE_COMPLEX(SERSIC_INDEX_);
    IDL_ENSURE_SCALAR(SERSIC_INDEX_);
    call[0] = SERSIC_INDEX_;
    SERSIC_INDEX_ = IDL_CvtDbl(1,call); /* May cause SERSIC_INDEX_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(SERSIC_INDEX_); /* Output cannot be expression */
    IDL_StoreScalarZero(SERSIC_INDEX_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* REFF_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(REFF_);
    IDL_ENSURE_SIMPLE(REFF_);
    IDL_EXCLUDE_STRING(REFF_);
    IDL_EXCLUDE_COMPLEX(REFF_);
    IDL_ENSURE_SCALAR(REFF_);
    call[0] = REFF_;
    REFF_ = IDL_CvtDbl(1,call); /* May cause REFF_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(REFF_); /* Output cannot be expression */
    IDL_StoreScalarZero(REFF_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* AXIAL_RATIO_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(AXIAL_RATIO_);
    IDL_ENSURE_SIMPLE(AXIAL_RATIO_);
    IDL_EXCLUDE_STRING(AXIAL_RATIO_);
    IDL_EXCLUDE_COMPLEX(AXIAL_RATIO_);
    IDL_ENSURE_SCALAR(AXIAL_RATIO_);
    call[0] = AXIAL_RATIO_;
    AXIAL_RATIO_ = IDL_CvtDbl(1,call); /* May cause AXIAL_RATIO_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(AXIAL_RATIO_); /* Output cannot be expression */
    IDL_StoreScalarZero(AXIAL_RATIO_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* BULGETODISK_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(BULGETODISK_);
    IDL_ENSURE_SIMPLE(BULGETODISK_);
    IDL_EXCLUDE_STRING(BULGETODISK_);
    IDL_EXCLUDE_COMPLEX(BULGETODISK_);
    IDL_ENSURE_SCALAR(BULGETODISK_);
    call[0] = BULGETODISK_;
    BULGETODISK_ = IDL_CvtDbl(1,call); /* May cause BULGETODISK_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(BULGETODISK_); /* Output cannot be expression */
    IDL_StoreScalarZero(BULGETODISK_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* INCLINATION_ANGLE_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(INCLINATION_ANGLE_);
    IDL_ENSURE_SIMPLE(INCLINATION_ANGLE_);
    IDL_EXCLUDE_STRING(INCLINATION_ANGLE_);
    IDL_EXCLUDE_COMPLEX(INCLINATION_ANGLE_);
    IDL_ENSURE_SCALAR(INCLINATION_ANGLE_);
    call[0] = INCLINATION_ANGLE_;
    INCLINATION_ANGLE_ = IDL_CvtDbl(1,call); /* May cause INCLINATION_ANGLE_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(INCLINATION_ANGLE_); /* Output cannot be expression */
    IDL_StoreScalarZero(INCLINATION_ANGLE_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* POSITION_ANGLE_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(POSITION_ANGLE_);
    IDL_ENSURE_SIMPLE(POSITION_ANGLE_);
    IDL_EXCLUDE_STRING(POSITION_ANGLE_);
    IDL_EXCLUDE_COMPLEX(POSITION_ANGLE_);
    IDL_ENSURE_SCALAR(POSITION_ANGLE_);
    call[0] = POSITION_ANGLE_;
    POSITION_ANGLE_ = IDL_CvtDbl(1,call); /* May cause POSITION_ANGLE_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(POSITION_ANGLE_); /* Output cannot be expression */
    IDL_StoreScalarZero(POSITION_ANGLE_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* PHASE_ANGLE_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(PHASE_ANGLE_);
    IDL_ENSURE_SIMPLE(PHASE_ANGLE_);
    IDL_EXCLUDE_STRING(PHASE_ANGLE_);
    IDL_EXCLUDE_COMPLEX(PHASE_ANGLE_);
    IDL_ENSURE_SCALAR(PHASE_ANGLE_);
    call[0] = PHASE_ANGLE_;
    PHASE_ANGLE_ = IDL_CvtDbl(1,call); /* May cause PHASE_ANGLE_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(PHASE_ANGLE_); /* Output cannot be expression */
    IDL_StoreScalarZero(PHASE_ANGLE_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DISTANCE_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DISTANCE_);
    IDL_ENSURE_SIMPLE(DISTANCE_);
    IDL_EXCLUDE_STRING(DISTANCE_);
    IDL_EXCLUDE_COMPLEX(DISTANCE_);
    IDL_ENSURE_SCALAR(DISTANCE_);
    call[0] = DISTANCE_;
    DISTANCE_ = IDL_CvtDbl(1,call); /* May cause DISTANCE_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DISTANCE_); /* Output cannot be expression */
    IDL_StoreScalarZero(DISTANCE_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* NXIM_ : INTEGER :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(NXIM_);
    IDL_ENSURE_SIMPLE(NXIM_);
    IDL_EXCLUDE_STRING(NXIM_);
    IDL_EXCLUDE_COMPLEX(NXIM_);
    IDL_ENSURE_SCALAR(NXIM_);
    call[0] = NXIM_;
    NXIM_ = IDL_CvtLng(1,call); /* May cause NXIM_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(NXIM_); /* Output cannot be expression */
    IDL_StoreScalarZero(NXIM_,IDL_TYP_LONG);  /* Not for in/out! */
  }

  in = 1;           /* NYIM_ : INTEGER :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(NYIM_);
    IDL_ENSURE_SIMPLE(NYIM_);
    IDL_EXCLUDE_STRING(NYIM_);
    IDL_EXCLUDE_COMPLEX(NYIM_);
    IDL_ENSURE_SCALAR(NYIM_);
    call[0] = NYIM_;
    NYIM_ = IDL_CvtLng(1,call); /* May cause NYIM_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(NYIM_); /* Output cannot be expression */
    IDL_StoreScalarZero(NYIM_,IDL_TYP_LONG);  /* Not for in/out! */
  }

  in = 1;           /* DXIM_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DXIM_);
    IDL_ENSURE_SIMPLE(DXIM_);
    IDL_EXCLUDE_STRING(DXIM_);
    IDL_EXCLUDE_COMPLEX(DXIM_);
    IDL_ENSURE_SCALAR(DXIM_);
    call[0] = DXIM_;
    DXIM_ = IDL_CvtDbl(1,call); /* May cause DXIM_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DXIM_); /* Output cannot be expression */
    IDL_StoreScalarZero(DXIM_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* DYIM_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(DYIM_);
    IDL_ENSURE_SIMPLE(DYIM_);
    IDL_EXCLUDE_STRING(DYIM_);
    IDL_EXCLUDE_COMPLEX(DYIM_);
    IDL_ENSURE_SCALAR(DYIM_);
    call[0] = DYIM_;
    DYIM_ = IDL_CvtDbl(1,call); /* May cause DYIM_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(DYIM_); /* Output cannot be expression */
    IDL_StoreScalarZero(DYIM_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 0;           /* IM_SCATT_ : REAL : (NXIM,NYIM) :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(IM_SCATT_);
    IDL_ENSURE_SIMPLE(IM_SCATT_);
    IDL_EXCLUDE_STRING(IM_SCATT_);
    IDL_EXCLUDE_COMPLEX(IM_SCATT_);
    IDL_ENSURE_ARRAY(IM_SCATT_);
    call[0] = IM_SCATT_;
    IM_SCATT_ = IDL_CvtFlt(1,call); /* May cause IM_SCATT_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(IM_SCATT_); /* Output cannot be expression */
    ndim = 2;
    dim[0] = NXIM_->value.l /*???*/;
    dim[1] = NYIM_->value.l /*???*/;
    IDL_StoreScalarZero(IM_SCATT_,IDL_TYP_FLOAT);  /* Free resources */
    IDL_MakeTempArray(IDL_TYP_FLOAT,ndim,dim,IDL_ARR_INI_ZERO,&tmp);
    IDL_VarCopy(tmp,IM_SCATT_);
  }

  in = 0;           /* IM_DIREC_ : REAL : (NXIM,NYIM) :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(IM_DIREC_);
    IDL_ENSURE_SIMPLE(IM_DIREC_);
    IDL_EXCLUDE_STRING(IM_DIREC_);
    IDL_EXCLUDE_COMPLEX(IM_DIREC_);
    IDL_ENSURE_ARRAY(IM_DIREC_);
    call[0] = IM_DIREC_;
    IM_DIREC_ = IDL_CvtFlt(1,call); /* May cause IM_DIREC_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(IM_DIREC_); /* Output cannot be expression */
    ndim = 2;
    dim[0] = NXIM_->value.l /*???*/;
    dim[1] = NYIM_->value.l /*???*/;
    IDL_StoreScalarZero(IM_DIREC_,IDL_TYP_FLOAT);  /* Free resources */
    IDL_MakeTempArray(IDL_TYP_FLOAT,ndim,dim,IDL_ARR_INI_ZERO,&tmp);
    IDL_VarCopy(tmp,IM_DIREC_);
  }

  in = 0;           /* IM_SCATT_SIG_ : REAL : (NXIM,NYIM) :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(IM_SCATT_SIG_);
    IDL_ENSURE_SIMPLE(IM_SCATT_SIG_);
    IDL_EXCLUDE_STRING(IM_SCATT_SIG_);
    IDL_EXCLUDE_COMPLEX(IM_SCATT_SIG_);
    IDL_ENSURE_ARRAY(IM_SCATT_SIG_);
    call[0] = IM_SCATT_SIG_;
    IM_SCATT_SIG_ = IDL_CvtFlt(1,call); /* May cause IM_SCATT_SIG_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(IM_SCATT_SIG_); /* Output cannot be expression */
    ndim = 2;
    dim[0] = NXIM_->value.l /*???*/;
    dim[1] = NYIM_->value.l /*???*/;
    IDL_StoreScalarZero(IM_SCATT_SIG_,IDL_TYP_FLOAT);  /* Free resources */
    IDL_MakeTempArray(IDL_TYP_FLOAT,ndim,dim,IDL_ARR_INI_ZERO,&tmp);
    IDL_VarCopy(tmp,IM_SCATT_SIG_);
  }

  in = 0;           /* IM_DIREC_SIG_ : REAL : (NXIM,NYIM) :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(IM_DIREC_SIG_);
    IDL_ENSURE_SIMPLE(IM_DIREC_SIG_);
    IDL_EXCLUDE_STRING(IM_DIREC_SIG_);
    IDL_EXCLUDE_COMPLEX(IM_DIREC_SIG_);
    IDL_ENSURE_ARRAY(IM_DIREC_SIG_);
    call[0] = IM_DIREC_SIG_;
    IM_DIREC_SIG_ = IDL_CvtFlt(1,call); /* May cause IM_DIREC_SIG_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(IM_DIREC_SIG_); /* Output cannot be expression */
    ndim = 2;
    dim[0] = NXIM_->value.l /*???*/;
    dim[1] = NYIM_->value.l /*???*/;
    IDL_StoreScalarZero(IM_DIREC_SIG_,IDL_TYP_FLOAT);  /* Free resources */
    IDL_MakeTempArray(IDL_TYP_FLOAT,ndim,dim,IDL_ARR_INI_ZERO,&tmp);
    IDL_VarCopy(tmp,IM_DIREC_SIG_);
  }

  in = 1;           /* PSF_FILE_ : CHARACTER*(128) :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(PSF_FILE_);
    IDL_ENSURE_SIMPLE(PSF_FILE_);
    IDL_EXCLUDE_COMPLEX(PSF_FILE_);
    IDL_ENSURE_SCALAR(PSF_FILE_);
    call[0] = PSF_FILE_;
    PSF_FILE_ = IDL_CvtByte(1,call); /* May cause PSF_FILE_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(PSF_FILE_); /* Output cannot be expression */
    IDL_StoreScalarZero(PSF_FILE_,IDL_TYP_BYTE);  /* Not for in/out! */
  }

  in = 1;           /* LEFT_RIGHT_FOLD_ : INTEGER :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(LEFT_RIGHT_FOLD_);
    IDL_ENSURE_SIMPLE(LEFT_RIGHT_FOLD_);
    IDL_EXCLUDE_STRING(LEFT_RIGHT_FOLD_);
    IDL_EXCLUDE_COMPLEX(LEFT_RIGHT_FOLD_);
    IDL_ENSURE_SCALAR(LEFT_RIGHT_FOLD_);
    call[0] = LEFT_RIGHT_FOLD_;
    LEFT_RIGHT_FOLD_ = IDL_CvtLng(1,call); /* May cause LEFT_RIGHT_FOLD_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(LEFT_RIGHT_FOLD_); /* Output cannot be expression */
    IDL_StoreScalarZero(LEFT_RIGHT_FOLD_,IDL_TYP_LONG);  /* Not for in/out! */
  }

     galaxy_idl_(GETVARADDR(NO_PHOTONS_),
                 GETVARADDR(NPRINT_),GETVARADDR(HGG_),
                 GETVARADDR(ALBEDO_),GETVARADDR(LUMINOSITY_),
                 GETVARADDR(DUST1_),GETVARADDR(TAUFACE1_),
                 GETVARADDR(DUST_RSCALE1_),
                 GETVARADDR(DUST_ZSCALE1_),
                 GETVARADDR(DUST_RMAX1_),
                 GETVARADDR(DUST_ZMAX1_),
                 GETVARADDR(DUST2_),GETVARADDR(TAUFACE2_),
                 GETVARADDR(DUST_RSCALE2_),
                 GETVARADDR(DUST_ZSCALE2_),
                 GETVARADDR(DUST_RMAX2_),
                 GETVARADDR(DUST_ZMAX2_),
                 GETVARADDR(RMAX_),GETVARADDR(PMAX_),
                 GETVARADDR(ZMAX_),GETVARADDR(NR_),
                 GETVARADDR(NP_),GETVARADDR(NZ_),
                 GETVARADDR(DISK_NAME_),GETVARADDR(DISK_RSCALE_),
                 GETVARADDR(DISK_ZSCALE_),
                 GETVARADDR(DISK_RMAX_),GETVARADDR(DISK_ZMAX_),
                 GETVARADDR(BULGE_NAME_),
                 GETVARADDR(SERSIC_INDEX_),
                 GETVARADDR(REFF_),GETVARADDR(AXIAL_RATIO_),
                 GETVARADDR(BULGETODISK_),
                 GETVARADDR(INCLINATION_ANGLE_),
                 GETVARADDR(POSITION_ANGLE_),
                 GETVARADDR(PHASE_ANGLE_),
                 GETVARADDR(DISTANCE_),GETVARADDR(NXIM_),
                 GETVARADDR(NYIM_),GETVARADDR(DXIM_),
                 GETVARADDR(DYIM_),GETVARADDR(IM_SCATT_),
                 GETVARADDR(IM_DIREC_),GETVARADDR(IM_SCATT_SIG_),
                 GETVARADDR(IM_DIREC_SIG_),
                 GETVARADDR(PSF_FILE_),GETVARADDR(LEFT_RIGHT_FOLD_));
  i=0;
  if (NO_PHOTONS_!=argv[i++]) IDL_DELTMP(NO_PHOTONS_);
  if (NPRINT_!=argv[i++]) IDL_DELTMP(NPRINT_);
  if (HGG_!=argv[i++]) IDL_DELTMP(HGG_);
  if (ALBEDO_!=argv[i++]) IDL_DELTMP(ALBEDO_);
  if (LUMINOSITY_!=argv[i++]) IDL_DELTMP(LUMINOSITY_);
  if (DUST1_!=argv[i++]) IDL_DELTMP(DUST1_);
  if (TAUFACE1_!=argv[i++]) IDL_DELTMP(TAUFACE1_);
  if (DUST_RSCALE1_!=argv[i++]) IDL_DELTMP(DUST_RSCALE1_);
  if (DUST_ZSCALE1_!=argv[i++]) IDL_DELTMP(DUST_ZSCALE1_);
  if (DUST_RMAX1_!=argv[i++]) IDL_DELTMP(DUST_RMAX1_);
  if (DUST_ZMAX1_!=argv[i++]) IDL_DELTMP(DUST_ZMAX1_);
  if (DUST2_!=argv[i++]) IDL_DELTMP(DUST2_);
  if (TAUFACE2_!=argv[i++]) IDL_DELTMP(TAUFACE2_);
  if (DUST_RSCALE2_!=argv[i++]) IDL_DELTMP(DUST_RSCALE2_);
  if (DUST_ZSCALE2_!=argv[i++]) IDL_DELTMP(DUST_ZSCALE2_);
  if (DUST_RMAX2_!=argv[i++]) IDL_DELTMP(DUST_RMAX2_);
  if (DUST_ZMAX2_!=argv[i++]) IDL_DELTMP(DUST_ZMAX2_);
  if (RMAX_!=argv[i++]) IDL_DELTMP(RMAX_);
  if (PMAX_!=argv[i++]) IDL_DELTMP(PMAX_);
  if (ZMAX_!=argv[i++]) IDL_DELTMP(ZMAX_);
  if (NR_!=argv[i++]) IDL_DELTMP(NR_);
  if (NP_!=argv[i++]) IDL_DELTMP(NP_);
  if (NZ_!=argv[i++]) IDL_DELTMP(NZ_);
  if (DISK_NAME_!=argv[i++]) IDL_DELTMP(DISK_NAME_);
  if (DISK_RSCALE_!=argv[i++]) IDL_DELTMP(DISK_RSCALE_);
  if (DISK_ZSCALE_!=argv[i++]) IDL_DELTMP(DISK_ZSCALE_);
  if (DISK_RMAX_!=argv[i++]) IDL_DELTMP(DISK_RMAX_);
  if (DISK_ZMAX_!=argv[i++]) IDL_DELTMP(DISK_ZMAX_);
  if (BULGE_NAME_!=argv[i++]) IDL_DELTMP(BULGE_NAME_);
  if (SERSIC_INDEX_!=argv[i++]) IDL_DELTMP(SERSIC_INDEX_);
  if (REFF_!=argv[i++]) IDL_DELTMP(REFF_);
  if (AXIAL_RATIO_!=argv[i++]) IDL_DELTMP(AXIAL_RATIO_);
  if (BULGETODISK_!=argv[i++]) IDL_DELTMP(BULGETODISK_);
  if (INCLINATION_ANGLE_!=argv[i++]) IDL_DELTMP(INCLINATION_ANGLE_);
  if (POSITION_ANGLE_!=argv[i++]) IDL_DELTMP(POSITION_ANGLE_);
  if (PHASE_ANGLE_!=argv[i++]) IDL_DELTMP(PHASE_ANGLE_);
  if (DISTANCE_!=argv[i++]) IDL_DELTMP(DISTANCE_);
  if (NXIM_!=argv[i++]) IDL_DELTMP(NXIM_);
  if (NYIM_!=argv[i++]) IDL_DELTMP(NYIM_);
  if (DXIM_!=argv[i++]) IDL_DELTMP(DXIM_);
  if (DYIM_!=argv[i++]) IDL_DELTMP(DYIM_);
  if (IM_SCATT_!=argv[i++]) IDL_DELTMP(IM_SCATT_);
  if (IM_DIREC_!=argv[i++]) IDL_DELTMP(IM_DIREC_);
  if (IM_SCATT_SIG_!=argv[i++]) IDL_DELTMP(IM_SCATT_SIG_);
  if (IM_DIREC_SIG_!=argv[i++]) IDL_DELTMP(IM_DIREC_SIG_);
  if (PSF_FILE_!=argv[i++]) IDL_DELTMP(PSF_FILE_);
  if (LEFT_RIGHT_FOLD_!=argv[i++]) IDL_DELTMP(LEFT_RIGHT_FOLD_);
  return NULL_VPTR;
}

int IDL_Load(void)
{

#if IDL_SINCE(5,3)

  static IDL_SYSFUN_DEF2 proc_def[] = {
    {(IDL_FUN_RET) GALAXY_IDL,"GALAXY_IDL",47,47,0,(void*)0}
    };
  return IDL_SysRtnAdd(proc_def,IDL_FALSE,1);

#else  /* Before v 5.3 */

  static IDL_SYSFUN_DEF proc_def[] = {
    {(IDL_FUN_RET) GALAXY_IDL,"GALAXY_IDL",47,47}
    };
  return IDL_AddSystemRoutine(proc_def,IDL_FALSE,1);

#endif

}
