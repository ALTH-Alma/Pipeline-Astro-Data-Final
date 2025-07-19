#include "fargo3d.h"

void _CondInit(int o) {
  
  int i,j,k;
  real r, omega;
  
  real *rho  = Density->field_cpu;
  real *cs   = Energy->field_cpu;
  real *vphi = Vx->field_cpu;
  real *vr   = Vy->field_cpu;
  
  real rhog, rhod;
  real vk;
  real size;
  real eps;
  real ff;
  
  i = j = k = 0;
  
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
	
	r     = Ymed(j);
	omega = sqrt(G*MSTAR/r/r/r);                       //Keplerian frequency
	rhog  = SIGMA0*pow(r/R0,-SIGMASLOPE);
	rhog *= exp(-pow(r/R0/CUTOFF,4));
#ifdef SIZEDRAG
	size = AMAX;
	eps =  EPSILON;
	rhod = eps*rhog;
#else
	rhod  = rhog*EPSILON;                              //Dust surface density
#endif
	
	if (Fluidtype == GAS) {
	  rho[l]   = rhog;
	  vphi[l]  = omega*r;
	  vr[l]    = 0.0;
	  cs[l]    = ASPECTRATIO;
	}
	
	if (Fluidtype == DUST) {
	  rho[l]  = rhod;
	  vphi[l] = omega*r;
	  vr[l]   = 0.0;
	  cs[l]   = 0.0;
	}
	
	vphi[l] -= OMEGAFRAME*r;
	
      }
    }
  }
}

void CondInit() {
#ifdef SIZEDRAG
  RHO_SOLID *= MSTAR/(MFACTOR*MSTAR_CGS) / pow(R0/(RFACTOR*R0_CGS),3);
  AMIN *= R0/(RFACTOR*R0_CGS);
  AMAX *= R0/(RFACTOR*R0_CGS);
#endif
  
  int id_gas = 0;
  int feedback = YES;
  //We first create the gaseous fluid and store it in the array Fluids[]
  Fluids[id_gas] = CreateFluid("gas",GAS);

  //We now select the fluid
  SelectFluid(id_gas);

  //and fill its fields
  _CondInit(id_gas);


#ifdef SIZEDRAG
  int id_dust;
  char dust_name[MAXNAMELENGTH];
  FILE *stream;
  char fullname[MAXLINELENGTH];
  sprintf (fullname, "%s/dust.par", OUTPUTDIR);
  stream = fopen_prs (fullname,"w");
  for(id_dust=1;id_dust<NFLUIDS;id_dust++){
    sprintf(dust_name,"dust%d",id_dust);
    Fluids[id_dust] = CreateFluid(dust_name,DUST);
    SelectFluid(id_dust);
    INSPECT_INT(id_dust);
    _CondInit(id_dust);
    real a,rs;
    a = AMAX;
    fprintf(stream, "%d\t%.15g\n", id_dust, a/(R0/(RFACTOR*R0_CGS)));
    rs = RHO_SOLID;
    ColRate(1.0/(a*rs), 0, id_dust,YES);
    }
  fclose(stream);
#else
  //We repeat the process for the dust fluids
  char dust_name[MAXNAMELENGTH];
  int id_dust;

  for(id_dust = 1; id_dust<NFLUIDS; id_dust++) {
    sprintf(dust_name,"dust%d",id_dust); //We assign different names to the dust fluids

    Fluids[id_dust]  = CreateFluid(dust_name, DUST);
    SelectFluid(id_dust);
    _CondInit();

  }

  /*We now fill the collision matrix (Feedback from dust included)
   Note: ColRate() moves the collision matrix to the device.
   If feedback=NO, gas does not feel the drag force.*/
  
  ColRate(INVSTOKES1, id_gas, 1, feedback);
  ColRate(INVSTOKES2, id_gas, 2, feedback);
  ColRate(INVSTOKES3, id_gas, 3, feedback);
#endif

#ifdef GPU
  DevMemcpyH2D(Alpha_d,Alpha,sizeof(real)*NFLUIDS*NFLUIDS);
#endif
}
