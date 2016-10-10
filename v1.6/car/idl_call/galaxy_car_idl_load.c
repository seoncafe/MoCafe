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
extern  galaxy_car_idl_(void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*);

    
IDL_VPTR GALAXY_CAR_IDL(int argc, IDL_VPTR argv[])
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
  IDL_VPTR ISEED4_=argv[i++];  /* ISEED4_ : INTEGER :  :in? */
  IDL_VPTR DUST1_=argv[i++];  /* DUST1_ : CHARACTER*(100) :  :in? */
  IDL_VPTR TAUFACE1_=argv[i++];  /* TAUFACE1_ : REAL*8 :  :in? */
  IDL_VPTR DUST_RSCALE1_=argv[i++];  /* DUST_RSCALE1_ : REAL*8 :  :in? */
  IDL_VPTR DUST_ZSCALE1_=argv[i++];  /* DUST_ZSCALE1_ : REAL*8 :  :in? */
  IDL_VPTR DUST_RMAX1_=argv[i++];  /* DUST_RMAX1_ : REAL*8 :  :in? */
  IDL_VPTR DUST_ZMAX1_=argv[i++];  /* DUST_ZMAX1_ : REAL*8 :  :in? */
  IDL_VPTR DUST2_=argv[i++];  /* DUST2_ : CHARACTER*(100) :  :in? */
  IDL_VPTR TAUFACE2_=argv[i++];  /* TAUFACE2_ : REAL*8 :  :in? */
  IDL_VPTR DUST_RSCALE2_=argv[i++];  /* DUST_RSCALE2_ : REAL*8 :  :in? */
  IDL_VPTR DUST_ZSCALE2_=argv[i++];  /* DUST_ZSCALE2_ : REAL*8 :  :in? */
  IDL_VPTR DUST_RMAX2_=argv[i++];  /* DUST_RMAX2_ : REAL*8 :  :in? */
  IDL_VPTR DUST_ZMAX2_=argv[i++];  /* DUST_ZMAX2_ : REAL*8 :  :in? */
  IDL_VPTR XMAX_=argv[i++];  /* XMAX_ : REAL*8 :  :in? */
  IDL_VPTR YMAX_=argv[i++];  /* YMAX_ : REAL*8 :  :in? */
  IDL_VPTR ZMAX_=argv[i++];  /* ZMAX_ : REAL*8 :  :in? */
  IDL_VPTR NX_=argv[i++];  /* NX_ : INTEGER :  :in? */
  IDL_VPTR NY_=argv[i++];  /* NY_ : INTEGER :  :in? */
  IDL_VPTR NZ_=argv[i++];  /* NZ_ : INTEGER :  :in? */
  IDL_VPTR SOURCE_NAME_=argv[i++];  /* SOURCE_NAME_ : CHARACTER*(100) :  :in? */
  IDL_VPTR SOURCE_RSCALE_=argv[i++];  /* SOURCE_RSCALE_ : REAL*8 :  :in? */
  IDL_VPTR SOURCE_ZSCALE_=argv[i++];  /* SOURCE_ZSCALE_ : REAL*8 :  :in? */
  IDL_VPTR SOURCE_RMAX_=argv[i++];  /* SOURCE_RMAX_ : REAL*8 :  :in? */
  IDL_VPTR SOURCE_ZMAX_=argv[i++];  /* SOURCE_ZMAX_ : REAL*8 :  :in? */
  IDL_VPTR OBS_ANGLES_=argv[i++];  /* OBS_ANGLES_ : REAL*8 : (3) :in? */
  IDL_VPTR OBS_DISTANCE_=argv[i++];  /* OBS_DISTANCE_ : REAL*8 :  :in? */
  IDL_VPTR NXIM_=argv[i++];  /* NXIM_ : INTEGER :  :in? */
  IDL_VPTR NYIM_=argv[i++];  /* NYIM_ : INTEGER :  :in? */
  IDL_VPTR DXIM_=argv[i++];  /* DXIM_ : REAL*8 :  :in? */
  IDL_VPTR DYIM_=argv[i++];  /* DYIM_ : REAL*8 :  :in? */
  IDL_VPTR IM_SCATT_=argv[i++];  /* IM_SCATT_ : REAL : (NXIM,NYIM) :in? */
  IDL_VPTR IM_DIREC_=argv[i++];  /* IM_DIREC_ : REAL : (NXIM,NYIM) :in? */



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

  in = 1;           /* ISEED4_ : INTEGER :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(ISEED4_);
    IDL_ENSURE_SIMPLE(ISEED4_);
    IDL_EXCLUDE_STRING(ISEED4_);
    IDL_EXCLUDE_COMPLEX(ISEED4_);
    IDL_ENSURE_SCALAR(ISEED4_);
    call[0] = ISEED4_;
    ISEED4_ = IDL_CvtLng(1,call); /* May cause ISEED4_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(ISEED4_); /* Output cannot be expression */
    IDL_StoreScalarZero(ISEED4_,IDL_TYP_LONG);  /* Not for in/out! */
  }

  in = 1;           /* DUST1_ : CHARACTER*(100) :  :  */
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

  in = 1;           /* DUST2_ : CHARACTER*(100) :  :  */
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

  in = 1;           /* XMAX_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(XMAX_);
    IDL_ENSURE_SIMPLE(XMAX_);
    IDL_EXCLUDE_STRING(XMAX_);
    IDL_EXCLUDE_COMPLEX(XMAX_);
    IDL_ENSURE_SCALAR(XMAX_);
    call[0] = XMAX_;
    XMAX_ = IDL_CvtDbl(1,call); /* May cause XMAX_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(XMAX_); /* Output cannot be expression */
    IDL_StoreScalarZero(XMAX_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* YMAX_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(YMAX_);
    IDL_ENSURE_SIMPLE(YMAX_);
    IDL_EXCLUDE_STRING(YMAX_);
    IDL_EXCLUDE_COMPLEX(YMAX_);
    IDL_ENSURE_SCALAR(YMAX_);
    call[0] = YMAX_;
    YMAX_ = IDL_CvtDbl(1,call); /* May cause YMAX_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(YMAX_); /* Output cannot be expression */
    IDL_StoreScalarZero(YMAX_,IDL_TYP_DOUBLE);  /* Not for in/out! */
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

  in = 1;           /* NX_ : INTEGER :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(NX_);
    IDL_ENSURE_SIMPLE(NX_);
    IDL_EXCLUDE_STRING(NX_);
    IDL_EXCLUDE_COMPLEX(NX_);
    IDL_ENSURE_SCALAR(NX_);
    call[0] = NX_;
    NX_ = IDL_CvtLng(1,call); /* May cause NX_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(NX_); /* Output cannot be expression */
    IDL_StoreScalarZero(NX_,IDL_TYP_LONG);  /* Not for in/out! */
  }

  in = 1;           /* NY_ : INTEGER :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(NY_);
    IDL_ENSURE_SIMPLE(NY_);
    IDL_EXCLUDE_STRING(NY_);
    IDL_EXCLUDE_COMPLEX(NY_);
    IDL_ENSURE_SCALAR(NY_);
    call[0] = NY_;
    NY_ = IDL_CvtLng(1,call); /* May cause NY_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(NY_); /* Output cannot be expression */
    IDL_StoreScalarZero(NY_,IDL_TYP_LONG);  /* Not for in/out! */
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

  in = 1;           /* SOURCE_NAME_ : CHARACTER*(100) :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(SOURCE_NAME_);
    IDL_ENSURE_SIMPLE(SOURCE_NAME_);
    IDL_EXCLUDE_COMPLEX(SOURCE_NAME_);
    IDL_ENSURE_SCALAR(SOURCE_NAME_);
    call[0] = SOURCE_NAME_;
    SOURCE_NAME_ = IDL_CvtByte(1,call); /* May cause SOURCE_NAME_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(SOURCE_NAME_); /* Output cannot be expression */
    IDL_StoreScalarZero(SOURCE_NAME_,IDL_TYP_BYTE);  /* Not for in/out! */
  }

  in = 1;           /* SOURCE_RSCALE_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(SOURCE_RSCALE_);
    IDL_ENSURE_SIMPLE(SOURCE_RSCALE_);
    IDL_EXCLUDE_STRING(SOURCE_RSCALE_);
    IDL_EXCLUDE_COMPLEX(SOURCE_RSCALE_);
    IDL_ENSURE_SCALAR(SOURCE_RSCALE_);
    call[0] = SOURCE_RSCALE_;
    SOURCE_RSCALE_ = IDL_CvtDbl(1,call); /* May cause SOURCE_RSCALE_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(SOURCE_RSCALE_); /* Output cannot be expression */
    IDL_StoreScalarZero(SOURCE_RSCALE_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* SOURCE_ZSCALE_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(SOURCE_ZSCALE_);
    IDL_ENSURE_SIMPLE(SOURCE_ZSCALE_);
    IDL_EXCLUDE_STRING(SOURCE_ZSCALE_);
    IDL_EXCLUDE_COMPLEX(SOURCE_ZSCALE_);
    IDL_ENSURE_SCALAR(SOURCE_ZSCALE_);
    call[0] = SOURCE_ZSCALE_;
    SOURCE_ZSCALE_ = IDL_CvtDbl(1,call); /* May cause SOURCE_ZSCALE_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(SOURCE_ZSCALE_); /* Output cannot be expression */
    IDL_StoreScalarZero(SOURCE_ZSCALE_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* SOURCE_RMAX_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(SOURCE_RMAX_);
    IDL_ENSURE_SIMPLE(SOURCE_RMAX_);
    IDL_EXCLUDE_STRING(SOURCE_RMAX_);
    IDL_EXCLUDE_COMPLEX(SOURCE_RMAX_);
    IDL_ENSURE_SCALAR(SOURCE_RMAX_);
    call[0] = SOURCE_RMAX_;
    SOURCE_RMAX_ = IDL_CvtDbl(1,call); /* May cause SOURCE_RMAX_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(SOURCE_RMAX_); /* Output cannot be expression */
    IDL_StoreScalarZero(SOURCE_RMAX_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* SOURCE_ZMAX_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(SOURCE_ZMAX_);
    IDL_ENSURE_SIMPLE(SOURCE_ZMAX_);
    IDL_EXCLUDE_STRING(SOURCE_ZMAX_);
    IDL_EXCLUDE_COMPLEX(SOURCE_ZMAX_);
    IDL_ENSURE_SCALAR(SOURCE_ZMAX_);
    call[0] = SOURCE_ZMAX_;
    SOURCE_ZMAX_ = IDL_CvtDbl(1,call); /* May cause SOURCE_ZMAX_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(SOURCE_ZMAX_); /* Output cannot be expression */
    IDL_StoreScalarZero(SOURCE_ZMAX_,IDL_TYP_DOUBLE);  /* Not for in/out! */
  }

  in = 1;           /* OBS_ANGLES_ : REAL*8 : (3) :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(OBS_ANGLES_);
    IDL_ENSURE_SIMPLE(OBS_ANGLES_);
    IDL_EXCLUDE_STRING(OBS_ANGLES_);
    IDL_EXCLUDE_COMPLEX(OBS_ANGLES_);
    IDL_ENSURE_ARRAY(OBS_ANGLES_);
    call[0] = OBS_ANGLES_;
    OBS_ANGLES_ = IDL_CvtDbl(1,call); /* May cause OBS_ANGLES_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(OBS_ANGLES_); /* Output cannot be expression */
    ndim = 1;
    dim[0] = 3 /* Should depend on input var? */;
    IDL_StoreScalarZero(OBS_ANGLES_,IDL_TYP_DOUBLE);  /* Free resources */
    IDL_MakeTempArray(IDL_TYP_DOUBLE,ndim,dim,IDL_ARR_INI_ZERO,&tmp);
    IDL_VarCopy(tmp,OBS_ANGLES_);
  }

  in = 1;           /* OBS_DISTANCE_ : REAL*8 :  :  */
  if (in) {
    IDL_EXCLUDE_UNDEF(OBS_DISTANCE_);
    IDL_ENSURE_SIMPLE(OBS_DISTANCE_);
    IDL_EXCLUDE_STRING(OBS_DISTANCE_);
    IDL_EXCLUDE_COMPLEX(OBS_DISTANCE_);
    IDL_ENSURE_SCALAR(OBS_DISTANCE_);
    call[0] = OBS_DISTANCE_;
    OBS_DISTANCE_ = IDL_CvtDbl(1,call); /* May cause OBS_DISTANCE_ to be tmp */
  } else { /* Output */ 
    IDL_EXCLUDE_EXPR(OBS_DISTANCE_); /* Output cannot be expression */
    IDL_StoreScalarZero(OBS_DISTANCE_,IDL_TYP_DOUBLE);  /* Not for in/out! */
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

     galaxy_car_idl_(GETVARADDR(NO_PHOTONS_),
                     GETVARADDR(NPRINT_),
                     GETVARADDR(HGG_),GETVARADDR(ALBEDO_),
                     GETVARADDR(LUMINOSITY_),
                     GETVARADDR(ISEED4_),
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
                     GETVARADDR(XMAX_),GETVARADDR(YMAX_),
                     GETVARADDR(ZMAX_),GETVARADDR(NX_),
                     GETVARADDR(NY_),GETVARADDR(NZ_),
                     GETVARADDR(SOURCE_NAME_),
                     GETVARADDR(SOURCE_RSCALE_),
                     GETVARADDR(SOURCE_ZSCALE_),
                     GETVARADDR(SOURCE_RMAX_),
                     GETVARADDR(SOURCE_ZMAX_),
                     GETVARADDR(OBS_ANGLES_),
                     GETVARADDR(OBS_DISTANCE_),
                     GETVARADDR(NXIM_),GETVARADDR(NYIM_),
                     GETVARADDR(DXIM_),GETVARADDR(DYIM_),
                     GETVARADDR(IM_SCATT_),
                     GETVARADDR(IM_DIREC_));
  i=0;
  if (NO_PHOTONS_!=argv[i++]) IDL_DELTMP(NO_PHOTONS_);
  if (NPRINT_!=argv[i++]) IDL_DELTMP(NPRINT_);
  if (HGG_!=argv[i++]) IDL_DELTMP(HGG_);
  if (ALBEDO_!=argv[i++]) IDL_DELTMP(ALBEDO_);
  if (LUMINOSITY_!=argv[i++]) IDL_DELTMP(LUMINOSITY_);
  if (ISEED4_!=argv[i++]) IDL_DELTMP(ISEED4_);
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
  if (XMAX_!=argv[i++]) IDL_DELTMP(XMAX_);
  if (YMAX_!=argv[i++]) IDL_DELTMP(YMAX_);
  if (ZMAX_!=argv[i++]) IDL_DELTMP(ZMAX_);
  if (NX_!=argv[i++]) IDL_DELTMP(NX_);
  if (NY_!=argv[i++]) IDL_DELTMP(NY_);
  if (NZ_!=argv[i++]) IDL_DELTMP(NZ_);
  if (SOURCE_NAME_!=argv[i++]) IDL_DELTMP(SOURCE_NAME_);
  if (SOURCE_RSCALE_!=argv[i++]) IDL_DELTMP(SOURCE_RSCALE_);
  if (SOURCE_ZSCALE_!=argv[i++]) IDL_DELTMP(SOURCE_ZSCALE_);
  if (SOURCE_RMAX_!=argv[i++]) IDL_DELTMP(SOURCE_RMAX_);
  if (SOURCE_ZMAX_!=argv[i++]) IDL_DELTMP(SOURCE_ZMAX_);
  if (OBS_ANGLES_!=argv[i++]) IDL_DELTMP(OBS_ANGLES_);
  if (OBS_DISTANCE_!=argv[i++]) IDL_DELTMP(OBS_DISTANCE_);
  if (NXIM_!=argv[i++]) IDL_DELTMP(NXIM_);
  if (NYIM_!=argv[i++]) IDL_DELTMP(NYIM_);
  if (DXIM_!=argv[i++]) IDL_DELTMP(DXIM_);
  if (DYIM_!=argv[i++]) IDL_DELTMP(DYIM_);
  if (IM_SCATT_!=argv[i++]) IDL_DELTMP(IM_SCATT_);
  if (IM_DIREC_!=argv[i++]) IDL_DELTMP(IM_DIREC_);
  return NULL_VPTR;
}

int IDL_Load(void)
{

#if IDL_SINCE(5,3)

  static IDL_SYSFUN_DEF2 proc_def[] = {
    {(IDL_FUN_RET) GALAXY_CAR_IDL,"GALAXY_CAR_IDL",37,37,0,(void*)0}
    };
  return IDL_SysRtnAdd(proc_def,IDL_FALSE,1);

#else  /* Before v 5.3 */

  static IDL_SYSFUN_DEF proc_def[] = {
    {(IDL_FUN_RET) GALAXY_CAR_IDL,"GALAXY_CAR_IDL",37,37}
    };
  return IDL_AddSystemRoutine(proc_def,IDL_FALSE,1);

#endif

}
