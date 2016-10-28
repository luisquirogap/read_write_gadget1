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
  float U;
  float rho;
  float Ne;
  float Nh;
  float h;
  float sfr;
  //float stellar_age;
  float ecr;
} gas_properties;

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
//GLOBAL VARIABLES
////////////////////////////////////////////////////////////////////////
particulas *particles;
gas_properties *gaspro;
io_header Header;

int N_part_total,N_part[6],N_min,N_max;
int n0,n1,n2,n3,n4,n5;
int returnRead;

double Mass_total,Mass_tot[6];


///////////////////////////////////////////////////////////////////////////////////////
//ROUTINES
///////////////////////////////////////////////////////////////////////////////////////
FILE *fileOpen(char filename[],char mode[]);
int gsl_int_int_sort(int dimension,int *fvector,int *ivector);

int main(int argc,char *argv[])
{
  
  int i,j,type;
  int nPartWithMass;
  FILE *fHeader,*fdata;
  char tmp[1000],*filename,outfile[1000];
  int dummy;
  double totalMasses[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};  
  
  ////////////////////////////////////////////////////////////////////////
  //READ DATA FILE
  ////////////////////////////////////////////////////////////////////////
  filename = argv[1];
  fdata = fileOpen(filename,"r");
  
  ////////////////////////////////////////////////////////////////////////
  //READ HEADER
  ////////////////////////////////////////////////////////////////////////
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  returnRead = fread(&Header,sizeof(io_header),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  ////////////////////////////////////////////////////////////////////////
  //SAVE HEADER
  ////////////////////////////////////////////////////////////////////////
  sprintf(tmp,"%s_Header.dat",filename);
  fHeader=fopen(tmp,"w");
  
  fprintf(fHeader,"\n");
  fprintf(fHeader,"Read header from %s\n",filename);
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
  
  
  ////////////////////////////////////////////////////////////////////////
  //TOTAL NUMBER OF PARTICLES
  ////////////////////////////////////////////////////////////////////////
  N_part_total = 0;
  nPartWithMass = 0;
  printf("Reading snapshot with:\n");
  for(i=0; i<6; i++)
    {
      N_part[i] = Header.npartTotal[i];
      N_part_total += N_part[i];
      printf("%d particles of type %d\n",N_part[i],i);
      if( Header.mass[i]>0.0 )
	nPartWithMass = nPartWithMass +0;
      else
	nPartWithMass = nPartWithMass + Header.npartTotal[i]; 
    }
  printf("%d particles in the snapshot\n",N_part_total);
  
  
  ////////////////////////////////////////////////////////////////////////
  //ALLOCATE AND READ
  ////////////////////////////////////////////////////////////////////////
  
  particles = (particulas *)malloc((size_t)N_part_total*sizeof(particulas));
  if(particles == NULL){
    printf("Allocation of particles failed\n");
    exit(0);
  }
  
  if( N_part[0] > 0 )
    {
      gaspro = (gas_properties *)malloc((size_t) N_part[0]*sizeof(gas_properties));
      if(gaspro == NULL){
	printf("Allocation of gaspro failed\n");
	exit(0);
      }
    }
  
  fprintf(fHeader,"\nRead Blocks:\n");
  
  //*************************************************************************
  //POSITIONS
  //*************************************************************************
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  for(i=0; i<N_part_total; i++)
    returnRead = fread(&particles[i].pos[0],sizeof(float),3,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  //*************************************************************************
  //VELOCITIES
  //*************************************************************************
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  for(i=0; i<N_part_total; i++)
    returnRead = fread(&particles[i].vel[0],sizeof(float),3,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  //*************************************************************************
  //IDS
  //*************************************************************************
  returnRead = fread(&dummy,sizeof(dummy),1,fdata); 
#ifdef LONGIDS
  for(i=0; i<N_part_total; i++)
    returnRead = fread(&particles[i].Id,sizeof(unsigned long long),1,fdata);
#else
  int indexmin=0, indexmax=0, counter=0;
  for(i=0;i<4;i++)
    indexmin = indexmin + N_part[i];
  indexmax = indexmin + N_part[4];
  for(i=0; i<N_part_total; i++)
    {
      returnRead = fread(&particles[i].id,sizeof(unsigned int),1,fdata);
      /*      //printf("%d %d %u\n",counter++,i,particles[i].id);
      
      //if( (i>=0) && (i<N_part[0]) )
      // printf("Gas %d %d %u\n",counter++,i,particles[i].id);
      if( (i>=indexmin) && (i<indexmax) )
	{
	  printf("SN %d %d %u\n",counter++,i,particles[i].id);
	  }*/
    }
#endif
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  //*************************************************************************
  //MASSES
  //*************************************************************************
  if( nPartWithMass>0  )  
    returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  N_min=N_max=0;
  
  for(j=0;j<=5;j++)
    {
      N_max=N_max+Header.npartTotal[j];
      if( Header.npartTotal[j]>0 )
	{
	  if( Header.mass[j]>0.0 )
	    {
	      for(i=N_min;i<N_max;i++)
		particles[i].mass = Header.mass[j];
	    }
	  else
	    {
	      for(i=N_min;i<N_max;i++)
		returnRead = fread(&particles[i].mass,sizeof(float),1,fdata); 
	    }
	  N_min=N_max;
	}
    }
  
  if( nPartWithMass>0  )  
    returnRead = fread(&dummy,sizeof(dummy),1,fdata);
     
  
  //int N_metal = N_part_total - Header.npartTotal[1];  
 
  //*****************************************************
  //Reading additional properties
  //*****************************************************

  if(Header.npartTotal[0]!=0)
    {
      
      //Read internal energy for particles of gas
      //===============================
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].U,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      //===============================
      
      //Read density for particles of gas
      //===============================
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].rho,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      //===============================
      
#ifdef COOLING
      //Read electron density for particles of gas
      //==============================
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].Ne,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      //==============================
      
      //Read neutral hydrigen density for particles of gas
      //==============================
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].Nh,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      //==============================
#endif
      
      //Read SPH smoothing length for particles of gas
      //===============================
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].h,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      //===============================
      
      
     	  //Read star formation rate for particles of gas
	  //=============================
#ifdef SFR
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].sfr,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif
		  //=============================
		  /*
		  if(Header.flag_sfr==1)
		    {
		      
		      //Read formation time of star for particles of gas
		      //=============================
	  if(stellarage==1)
	    {
	      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
	      for(i=0;i<Header.npartTotal[0];i++)
			  returnRead = fread(&gaspro[i].stellar_age,sizeof(float),1,fdata);
		      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
	    }
	  //=============================
	}
      //==============================
      */
    }
  /*
  //Read metallicity for gas and stars
  //===============================
  if((Header.flag_sfr==1) && (metallicity==1))
    {
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<N_metal;i++)
	returnRead = fread(&Z[i],sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
    }
  //===============================
  */

  //Read gravitational potential for all particles
  //===============================
#ifdef OUTPUTPOTENTIAL
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  for(i=0;i<N_part_total;i++)
    returnRead = fread(&particles[i].pot,sizeof(float),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif
  //===============================
  
  //Read acceleration for all particles
  //===============================
#ifdef OUTPUTACCELERATION
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  for(i=0;i<N_part_total;i++)
    returnRead = fread(&particles[i].acce,sizeof(float),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif
  //===============================
  
  //Read rate of change of entropic function for gas 
  //===============================
  if( N_part[0] >0 )
    {
#ifdef OUTPUTCHANGEOFENTROPY
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].ecr,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif
    }
  //===============================
  
  //Read timestep for all particles
  //=============================== 
#ifdef OUTPUTTIMESTEP
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  for(i=0;i<N_part_total;i++)
    returnRead = fread(&particles[i].timestep,sizeof(float),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif
  //====================================================
  
  fclose(fHeader);  
  
  //######################################
  //Sorting particles by id
  //######################################
  
  particulas  *aux;
  aux = (particulas *)malloc((size_t)N_part_total*sizeof(particulas));
  if(particles == NULL){
    printf("Allocation of aux failed\n");
    exit(0);
  }

  size_t  *p;
  p = (size_t *)malloc((size_t)N_part_total*sizeof(size_t));
  if(p == NULL){
    printf("Allocation of p failed\n");
    exit(0);
  }

#ifdef LONGIDS
  unsigned long long *ID;
  ID = (unsigned long long *)malloc((size_t)N_part_total*sizeof(unsigned long long));
  if(ID == NULL){
    printf("Allocation of ID failed\n");
    exit(0);
  }
#else
  unsigned int *ID;
  ID = (unsigned int *)malloc((size_t)N_part_total*sizeof(unsigned int));
  if(ID == NULL){
    printf("Allocation of ID failed\n");
    exit(0);
  }
#endif
  
  for( i=0; i<N_part_total; i++)
    {  
      ID[i] = particles[i].id;  
      aux[i] = particles[i];
    }
  
   gsl_sort_uint_index(p,ID,1,(size_t) N_part_total); 		
   
   for( i=0; i<N_part_total; i++)
    particles[i] = aux[p[i]];   
   
   
   //int counter=0;
   //for( i=0; i<N_part_total; i++)
   //{
   //  printf("%d %u %lu %u\n",i,ID[i],p[i],particles[i].id);
   //  if(particles[i].id>=400000)
   //	 counter++;
   //}
   //printf("%d Estrellas nuevas\n",counter);
   
   
  if( N_part[0] > 0 )
    {
      gas_properties *auxgas;
      auxgas = (gas_properties *)malloc((size_t) N_part[0]*sizeof(gas_properties));
      if(auxgas == NULL){
	printf("Allocation of auxgas failed\n");
	exit(0);
      }
      
      for( i=0; i<N_part[0]; i++)
	auxgas[i] = gaspro[i];

	
      //printf("%d \t",i);      

      for( i=0; i<N_part[0]; i++)
	{
	  // printf("%d %lu %.8f\n",i,p[i],particles[i].mass);
	  gaspro[i] = auxgas[p[i]];
	}      
      
      free(auxgas);
    }
  
  free(ID);
  free(aux);
  free(p);
  
  n0 = Header.npartTotal[0];
  n1 = Header.npartTotal[1];
  n2 = Header.npartTotal[2];
  n3 = Header.npartTotal[3];
  n4 = Header.npartTotal[4];
  n5 = Header.npartTotal[5];
  
  for( i=0; i<N_part_total; i++)
    {
      if(i < n0) 
	particles[i].type = 0;
      
      if( (i >= n0) && ( i < (n0+n1)) ) 
	particles[i].type = 1;
      
      if( (i >= (n0+n1)) && (i < (n0+n1+n2)) ) 
	particles[i].type = 2;
      
      if( (i >= (n0+n1+n2) ) && (i < (n0+n1+n2+n3)) ) 
	particles[i].type = 3;
      
      if( (i >= (n0+n1+n2+n3)) && (i < (n0+n1+n2+n3+n4)) ) 
	particles[i].type = 4;
      
      if((i >= (n0+n1+n2+n3+n4)) && (i < (n0+n1+n2+n3+n4+n5)) ) 
	particles[i].type = 5;
    }   

 
  
  //****************************************************
  // Writing data for all particles in ASCII 
  //**************************************************** 

  N_min = 0;
  totalMasses[7] = 0.0;
  
  for( type=0; type<6; type++)
    {

      N_max = N_min + N_part[type];	 

      if( N_part[type] > 0 )
	{
	  sprintf(outfile,"%s.%d",filename,type);
	  FILE *outfiles;
	  outfiles = fopen(outfile,"w");
	  if(outfiles==NULL) printf("No se pudo abrir %s\n",outfile);  
	  
	  for(i=N_min;i<N_max;i++)
	    {

	      totalMasses[type] = totalMasses[type] + particles[i].mass;  
	      
#ifdef LONGIDS
	      fprintf(outfiles,"%lu",particles[i].id);
#else
	      fprintf(outfiles,"%u",particles[i].id);
#endif
	      
	      fprintf(outfiles," %f %f %f %f %f %f %f",
		      particles[i].pos[0],particles[i].pos[1],particles[i].pos[2],
		      particles[i].vel[0],particles[i].vel[1],particles[i].vel[2],
		      particles[i].mass);
	      
	      if( type == 0 &&  N_part[0] >0 )
		{
		  fprintf(outfiles," %f %f",gaspro[i].U,gaspro[i].rho);
#ifdef COOLING
		  fprintf(outfiles," %f %f",gaspro[i].Ne,gaspro[i].Nh);
#endif
		  fprintf(outfiles," %f",gaspro[i].h);
#ifdef SFR
		  fprintf(outfiles," %f",gaspro[i].sfr);
#endif
		}
	      
#ifdef OUTPUTPOTENTIAL
	      fprintf(outfiles," %f",particles[i].pot);
#endif
	      
#ifdef OUTPUTACCELERATION
	      fprintf(outfiles," %f",particles[i].acce);
#endif
	      
	
#ifdef OUTPUTCHANGEOFENTROP
	      if( type == 0 && N_part[0] >0 )
		{
		  fprintf(outfiles," %f",gaspro[i].ecr);
		}
#endif 
	
	      
#ifdef OUTPUTTIMESTEP
	      fprintf(outfiles," %f",particles[i].timestep);
#endif
	      
	      fprintf(outfiles,"\n");	
	      
	    }
	  fclose(outfiles);
	}

      N_min = N_max;

      printf("Total mass type %d = %.8e\n",type,totalMasses[type]);
      totalMasses[7] = totalMasses[7] + totalMasses[type];
    }

  printf("Galaxy total mass = %.8e\n",totalMasses[7]);  

  printf("Done.\n\n");

  free(particles);
  free(gaspro);  
  
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

