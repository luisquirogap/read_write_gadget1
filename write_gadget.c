////////////////////////////////////////////////////////////////////////
//HEADERS
////////////////////////////////////////////////////////////////////////
#include<stdio.h>
#include<stdlib.h>  
#include<string.h>  
#include<malloc.h>
#include<math.h>
#include<time.h>

//CABECERAS DE GSL                                                                  
#include<gsl/gsl_errno.h>
#include<gsl/gsl_sort.h>
#include <gsl/gsl_permutation.h>
#include<gsl/gsl_sort_uint.h>

////////////////////
// MACROS 
///////////////////

#define X 0
#define Y 1
#define Z 2

////////////////////////////////////////////////////////////////////////
//DATA STRUCTURES
////////////////////////////////////////////////////////////////////////
typedef struct 
{
  int type;
#ifdef LONGIDS
  unsigned long long id;
#else
  unsigned int id;
#endif
  float pos[3];
  float vel[3];
  float mass;
  //float Z;
  float pot;
  float acce;
  float timestep;
}particulas;

typedef struct
{
  int Npart[6];
  double mass[6];
  double time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned int npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];
  int  flag_entropy_instead_u;
  char     fill[256 - 6*4 - 6*8 - 2*8 - 2*4 - 6*4 - 2*4 - 4*8 - 2*4 - 6*4 - 1*4];  /* fills to 256 Bytes */
} io_header;

////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES
////////////////////////////////////////////////////////////////////////
particulas *particles;
io_header Header;

int N_part_total,N_part[6],N_min,N_max;
int n0,n1,n2,n3,n4,n5;
int returnRead;

double Mass_total,Mass_tot[6];
float *U;

///////////////////////////////////////////////////////////////////////////////////////
// ROUTINES
///////////////////////////////////////////////////////////////////////////////////////
FILE *fileOpen(char filename[],char mode[]);
int gsl_int_int_sort(int dimension,int *fvector,int *ivector);

int main(int argc,char *argv[])
{
  
  int i,imin=0,type;
  FILE *fParam,*fHeader,*fData,*fGadget;
  char *filename,outfile[500],outfile2[500];
  int dummy;
    
  
  ////////////////////////////////////////////////////////////////////////
  // READ DATA FILE
  ////////////////////////////////////////////////////////////////////////
  filename = argv[1];
  fParam = fileOpen(filename,"r");
  
  ////////////////////////////////////////////////////////////////////////
  // NAME OUTPUT FILE
  ////////////////////////////////////////////////////////////////////////
  returnRead = fscanf(fParam,"%s",outfile);
  fData = fileOpen(outfile,"r");
  printf("\nBuilding Gadget binary with format 1 in:\n%s\n\n",outfile);

  sprintf(outfile2,"%s.gad",outfile);
  fGadget = fopen(outfile2,"w");

  
  ////////////////////////////////////////////////////////////////////////
  // READ HEADER
  ////////////////////////////////////////////////////////////////////////
  returnRead = fscanf(fParam,"%d %d %d %d %d %d",
	 &Header.Npart[0],&Header.Npart[1],&Header.Npart[2],
	 &Header.Npart[3],&Header.Npart[4],&Header.Npart[5]);
  returnRead = fscanf(fParam,"%lf %lf %lf %lf %lf %lf",
	 &Header.mass[0],&Header.mass[1],&Header.mass[2],
	 &Header.mass[3],&Header.mass[4],&Header.mass[5]);
  returnRead = fscanf(fParam,"%lf",&Header.time);
  returnRead = fscanf(fParam,"%lf",&Header.redshift);
  returnRead = fscanf(fParam,"%d",&Header.flag_sfr);
  returnRead = fscanf(fParam,"%d",&Header.flag_feedback);
  returnRead = fscanf(fParam,"%u %u %u %u %u %u",
	 &Header.npartTotal[0],&Header.npartTotal[1],&Header.npartTotal[2],
	 &Header.npartTotal[3],&Header.npartTotal[4],&Header.npartTotal[5]);
  returnRead = fscanf(fParam,"%d",&Header.flag_cooling);
  returnRead = fscanf(fParam,"%d",&Header.num_files);
  returnRead = fscanf(fParam,"%lf",&Header.BoxSize);
  returnRead = fscanf(fParam,"%lf",&Header.Omega0);
  returnRead = fscanf(fParam,"%lf",&Header.OmegaLambda);
  returnRead = fscanf(fParam,"%lf",&Header.HubbleParam);
  returnRead = fscanf(fParam,"%d",&Header.flag_stellarage);
  returnRead = fscanf(fParam,"%d",&Header.flag_metals);
  returnRead = fscanf(fParam,"%u %u %u %u %u %u",
	 &Header.npartTotalHighWord[0],&Header.npartTotalHighWord[1],
	 &Header.npartTotalHighWord[2],&Header.npartTotalHighWord[3],
	 &Header.npartTotalHighWord[4],&Header.npartTotalHighWord[5]);
 returnRead = fscanf(fParam,"%d",&Header.flag_entropy_instead_u);

 fclose(fParam);
  
 ////////////////////////////////////////////////////////////////////////
 // SAVE HEADER
 ////////////////////////////////////////////////////////////////////////
 sprintf(outfile,"%s_Header_wrote.output",outfile);
 fHeader=fopen(outfile,"w");
 
 fprintf(fHeader,"\n");
 fprintf(fHeader,"Wrote Header for %s\n",filename);
 fprintf(fHeader,"\n");
 fprintf(fHeader,"int Npart[6] : \n");
 fprintf(fHeader,"Npart[0] = %d \n",Header.Npart[0]);
 fprintf(fHeader,"Npart[1] = %d \n",Header.Npart[1]);
 fprintf(fHeader,"Npart[2] = %d \n",Header.Npart[2]);
 fprintf(fHeader,"Npart[3] = %d \n",Header.Npart[3]);
 fprintf(fHeader,"Npart[4] = %d \n",Header.Npart[4]);
 fprintf(fHeader,"Npart[5] = %d \n",Header.Npart[5]);
 fprintf(fHeader,"double mass[6] : \n");
 fprintf(fHeader,"mass[0] = %.10e \n",Header.mass[0]);
 fprintf(fHeader,"mass[1] = %.10e \n",Header.mass[1]);
 fprintf(fHeader,"mass[2] = %.10e \n",Header.mass[2]);
 fprintf(fHeader,"mass[3] = %.10e \n",Header.mass[3]);
 fprintf(fHeader,"mass[4] = %.10e \n",Header.mass[4]);
 fprintf(fHeader,"mass[5] = %.10e \n",Header.mass[5]);
 fprintf(fHeader,"double time = %lf \n",Header.time);
 fprintf(fHeader,"double redshift = %lf \n",Header.redshift);
 fprintf(fHeader,"int flag_sfr = %d \n",Header.flag_sfr);
 fprintf(fHeader,"int flag_feedback = %d \n",Header.flag_feedback);
 fprintf(fHeader,"int npartTotal[6] : \n");
 fprintf(fHeader,"npartTotal[0] = %u \n",Header.npartTotal[0]);
 fprintf(fHeader,"npartTotal[1] = %u \n",Header.npartTotal[1]);
 fprintf(fHeader,"npartTotal[2] = %u \n",Header.npartTotal[2]);
 fprintf(fHeader,"npartTotal[3] = %u \n",Header.npartTotal[3]);
 fprintf(fHeader,"npartTotal[4] = %u \n",Header.npartTotal[4]);
 fprintf(fHeader,"npartTotal[5] = %u \n",Header.npartTotal[5]);
 fprintf(fHeader,"int flag_cooling = %d \n",Header.flag_cooling);
 fprintf(fHeader,"int num_files = %d \n",Header.num_files); 
 fprintf(fHeader,"double BoxSize = %lf \n",Header.BoxSize);
 fprintf(fHeader,"double Omega0 = %lf \n",Header.Omega0);
 fprintf(fHeader,"double OmegaLambda = %lf \n",Header.OmegaLambda);
 fprintf(fHeader,"double HubbleParam = %lf \n",Header.HubbleParam);
 fprintf(fHeader,"int flag_stellarage = %d \n",Header.flag_stellarage);
 fprintf(fHeader,"int flag_metals = %d \n",Header.flag_metals);
 fprintf(fHeader,"int npartTotalHighWord[6] : \n");
 fprintf(fHeader,"npartTotalHighWord[0] = %u \n",Header.npartTotalHighWord[0]);
 fprintf(fHeader,"npartTotalHighWord[1] = %u \n",Header.npartTotalHighWord[1]);
 fprintf(fHeader,"npartTotalHighWord[2] = %u \n",Header.npartTotalHighWord[2]);
 fprintf(fHeader,"npartTotalHighWord[3] = %u \n",Header.npartTotalHighWord[3]);
 fprintf(fHeader,"npartTotalHighWord[4] = %u \n",Header.npartTotalHighWord[4]);
 fprintf(fHeader,"npartTotalHighWord[5] = %u \n",Header.npartTotalHighWord[5]);
 fprintf(fHeader,"int flag_entropy_instead_u = %d \n",Header.flag_entropy_instead_u);
 fprintf(fHeader,"Rest size to 256 Bytes = %lu",sizeof(Header.fill));

 fclose(fHeader);   	       
 
 ////////////////////////////////////////////////////////////////////////
 // TOTAL NUMBER OF PARTICLES
 ////////////////////////////////////////////////////////////////////////
 N_part_total = 0;
 for(i=0; i<6; i++)
   {
     N_part[i] = Header.npartTotal[i];
     N_part_total += N_part[i];
   }
 
 ////////////////////////////////////////////////////////////////////////
 // ALLOCATE AND READ
 ////////////////////////////////////////////////////////////////////////
 
 particles = (particulas *)malloc((size_t)N_part_total*sizeof(particulas));
 if(particles == NULL){
   printf("Allocation of particles failed\n");
   exit(0);
 }
 
 U = (float *)malloc((size_t)Header.npartTotal[0]*sizeof(float));
 if(particles == NULL){
   printf("Allocation of U failed\n");
   exit(0);
 }
 
 ////////////////////////////////////////////////////////////
 // READING PARTICLES DATA
 ///////////////////////////////////////////////////////////
 
 //////////////////////////////////
 // READING GAS
 //////////////////////////////////
 if( Header.npartTotal[0] != 0)
   {
     for( i=0; i<Header.npartTotal[0]; i++)
       {
#ifdef LONGIDS
	 returnRead = fscanf(fData,"%llu %f %f %f %f %f %f %f %f",
		&particles[i].id,
		&particles[i].pos[X],&particles[i].pos[Y],&particles[i].pos[Z],
		&particles[i].vel[X],&particles[i].vel[Y],&particles[i].vel[Z],
		&particles[i].mass,
		&U[i]);
#else
	 returnRead = fscanf(fData,"%u %f %f %f %f %f %f %f %f",
		&particles[i].id,
		&particles[i].pos[X],&particles[i].pos[Y],&particles[i].pos[Z],
		&particles[i].vel[X],&particles[i].vel[Y],&particles[i].vel[Z],
		&particles[i].mass,
		&U[i]);
#endif
       }
     
     imin = Header.npartTotal[0];
   }
 
////////////////////////////////////////////////////////
// READING OTHER PARTICLES TYPE
////////////////////////////////////////////////////////
 for( i=imin; i<N_part_total; i++)
   {
#ifdef LONGIDS
     returnRead = fscanf(fData,"%llu %f %f %f %f %f %f %f",
	    &particles[i].id,
	    &particles[i].pos[X],&particles[i].pos[Y],&particles[i].pos[Z],
	    &particles[i].vel[X],&particles[i].vel[Y],&particles[i].vel[Z],
	    &particles[i].mass);
#else
     returnRead = fscanf(fData,"%u %f %f %f %f %f %f %f",
	    &particles[i].id,
	    &particles[i].pos[X],&particles[i].pos[Y],&particles[i].pos[Z],
	    &particles[i].vel[X],&particles[i].vel[Y],&particles[i].vel[Z],
	    &particles[i].mass);
#endif
   }
 
 fclose(fData);

 //////////////////////////////////
 // WRITING HEADER
 /////////////////////////////////
 dummy = sizeof(Header);
 fwrite(&dummy,sizeof(dummy),1,fGadget);
 fwrite(&Header,sizeof(Header),1,fGadget);
 fwrite(&dummy,sizeof(dummy),1,fGadget);

 ////////////////////////////////////////////////////////////////////////
 // WRITING POSITIONS
 ////////////////////////////////////////////////////////////////////////
 dummy = 3*N_part_total*sizeof(float);
 fwrite(&dummy,sizeof(dummy),1,fGadget); 
 for( i=0; i<N_part_total; i++ ) 
   fwrite(&particles[i].pos,sizeof(float),3,fGadget);
 fwrite(&dummy,sizeof(dummy),1,fGadget);

 ////////////////////////////////////////////////////////////////////////
 // WRITING VELOCITIES
 ////////////////////////////////////////////////////////////////////////
 dummy = 3*N_part_total*sizeof(float);
 fwrite(&dummy,sizeof(dummy),1,fGadget); 
 for( i=0; i<N_part_total; i++ ) 
   fwrite(&particles[i].vel,sizeof(float),3,fGadget);
 fwrite(&dummy,sizeof(dummy),1,fGadget);

 ////////////////////////////////////////////////////////////////////////
 // WRITING IDs
 ////////////////////////////////////////////////////////////////////////
#ifdef LONGIDS
 dummy = N_part_total*sizeof(unsigned long long);	
#else
 dummy = N_part_total*sizeof(unsigned int);		
#endif

 fwrite(&dummy,sizeof(dummy),1,fGadget);
#ifdef LONGIDS
 for(i=0; i<N_part_total; i++)
   fwrite(&particles[i].Id,sizeof(unsigned long long),1,fGadget);
#else
 for(i=0; i<N_part_total; i++)
   fwrite(&particles[i].id,sizeof(unsigned int),1,fGadget);
#endif
 fwrite(&dummy,sizeof(dummy),1,fGadget);
 
 ////////////////////////////////////////////////
 // WRITING MASSES
 ///////////////////////////////////////////////
 dummy = 0; 
 for( type=0; type<6; type++)
   if( Header.npartTotal[type]>0 )
     {
       if( Header.mass[type]>0.0 )
	 continue;
       else
	 dummy = dummy + Header.npartTotal[type]*sizeof(float);
     } 
 
 if(dummy>0)
   {
     fwrite(&dummy,sizeof(dummy),1,fGadget);
     
     N_min = N_max = 0;
     for( type=0; type<6; type++)
       {
	 N_max = N_max + Header.npartTotal[type];
	 if( (Header.npartTotal[type]>0) && (Header.mass[type]>0.0) )
	   continue;
	 else
	   {
	     for( i=N_min; i<N_max; i++)
	       {
		 fwrite(&particles[i].mass,sizeof(float),1,fGadget);
	       }
	     N_min = N_max;
	   }
       } 
     fwrite(&dummy,sizeof(dummy),1,fGadget);
   }

 ////////////////////////////////////////////////////////////
 // WRITING INTERNAL ENERGY FOR GAS
 ///////////////////////////////////////////////////////////
 if( Header.npartTotal[0]>0 )
   {
     dummy = Header.npartTotal[0]*sizeof(float);
     fwrite(&dummy,sizeof(dummy),1,fGadget); 
     for( i=0; i<Header.npartTotal[0]; i++ ) 
       fwrite(&U[i],sizeof(float),1,fGadget);
     fwrite(&dummy,sizeof(dummy),1,fGadget);
   }

 printf("Initial conditions in %s\n\n",outfile2);

 fclose(fGadget);   	      
  
 free(particles);
 free(U);
 
 return 0;
}

/*
  ======================================================================
  Open a file
  ======================================================================
*/
FILE *fileOpen(char filename[],char mode[])
{
  FILE *f;

  if( !(f=fopen(filename,mode)) ){
    fprintf(stderr,"Error openning '%s' for %s\n",filename,mode);
    exit(1);
  }
  
  return f;
}

