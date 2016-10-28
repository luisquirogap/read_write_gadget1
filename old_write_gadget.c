#include<stdio.h>
#include<stdlib.h>  
#include<malloc.h>

struct io_header
{
  int Npart[6];
  double mass[6];
  double time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
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
  char     fill[256 - 6*4 - 6*8 - 2*8 - 2*4 - 6*4 - 2*4 - 4*8 - 2*4 - 6*4 - 1*4];
} header;

struct nocolisional
{
  int id;
  float pos[3];
  float vel[3];
  float masa;
};

struct nocolisional *part[5];

struct colisional
{
  int id;
  float pos[3];
  float vel[3];
  float masa;
  float U;
  float rho;
  float Ne;
  float Neutral1H;
  float eps;
  float sfr;
};

struct colisional *gas;

int main(int argc,char *argv[])
{

  int i,j,k;
  int dummi,Npart;
  float fac_long,fac_masa,fac_vel;
  char *name;

  name=argv[2];

  FILE *outfile;
  outfile=fopen(name,"w");

  FILE *h;
  h=fopen("header.dat","r");

  fscanf(h,"%d %d %d %d %d %d",&header.Npart[0],&header.Npart[1],&header.Npart[2],
	 &header.Npart[3],&header.Npart[4],&header.Npart[5]);
  fscanf(h,"%lf %lf %lf %lf %lf %lf",&header.mass[0],&header.mass[1],&header.mass[2],
	 &header.mass[3],&header.mass[4],&header.mass[5]);
  fscanf(h,"%lf %lf %d %d",&header.time,&header.redshift,&header.flag_sfr,&header.flag_feedback);
  fscanf(h,"%d %d %d %d %d %d",&header.npartTotal[0],&header.npartTotal[1],&header.npartTotal[2],
	 &header.npartTotal[3],&header.npartTotal[4],&header.npartTotal[5]);
  fscanf(h,"%d %d %lf %lf %lf %lf",&header.flag_cooling,&header.num_files,
	 &header.BoxSize,&header.Omega0,&header.OmegaLambda,&header.HubbleParam);

  printf("%d %d %d %d %d %d\n",header.Npart[0],header.Npart[1],header.Npart[2],
	 header.Npart[3],header.Npart[4],header.Npart[5]);
  printf("%lf %lf %lf %lf %lf %lf\n",header.mass[0],header.mass[1],header.mass[2],
	 header.mass[3],header.mass[4],header.mass[5]);
  printf("%lf %lf %d %d\n",header.time,header.redshift,header.flag_sfr,header.flag_feedback);
  printf("%d %d %d %d %d %d\n",header.npartTotal[0],header.npartTotal[1],header.npartTotal[2],
	 header.npartTotal[3],header.npartTotal[4],header.npartTotal[5]);
  printf("%d %d %lf %lf %lf %lf\n",header.flag_cooling,header.num_files,
	 header.BoxSize,header.Omega0,header.OmegaLambda,header.HubbleParam);

  fclose(h);

  //  fac_long=0.286;
  // fac_masa=0.176;
  // fac_vel=0.00382;

  fac_long=1.0;
  fac_masa=1.0;
  fac_vel=1.0;


  printf("Factor longitud %f\n",fac_long);
  printf("Factor masa %f\n",fac_masa);
  printf("Factor velocidad %f\n",fac_vel);

  printf("%d %d %d %d\n",sizeof(int),sizeof(header.Npart[0]),sizeof(double),sizeof(header.mass[0]));

  Npart=header.npartTotal[0]+header.npartTotal[1]+header.npartTotal[2]+header.npartTotal[3]+header.npartTotal[4]+header.npartTotal[5];

  dummi=256;

  printf("%d %d %d\n",sizeof(header),dummi,Npart);
  //Write the header in the gadget format read from header.dat
  //==========================================================
  fwrite(&dummi,sizeof(dummi),1,outfile);
  fwrite(&header,sizeof(header),1,outfile);
  fwrite(&dummi,sizeof(dummi),1,outfile);
  //==========================================================

  //====================================
  //Note: All datas are read from infile
  //====================================

  //Charge all datas 
  //=======================================
  char names[100];
 

  for(k=0;k<6;k++)
    {	
    
      //if(header.npartTotal[k]>0)
      //{
	  sprintf(names,"%s.%d",argv[1],k);
	  printf("%s\n",names);
	  FILE *infile;	 
	  infile=fopen(names,"r");
	  //	}
      if(k==0)
	{
	  if(header.npartTotal[0]>0)
	    {
	      gas=malloc(header.npartTotal[0]*sizeof(struct colisional));
	      for(i=0;i<header.npartTotal[0];i++)
		{
		  fscanf(infile,"%d %f %f %f %f %f %f %f %f %f %f %f %f %f",
			 &gas[i].id,
			 &gas[i].pos[0],&gas[i].pos[1],&gas[i].pos[2],
			 &gas[i].vel[0],&gas[i].vel[1],&gas[i].vel[2],
			 &gas[i].masa,
			 &gas[i].U,
			 &gas[i].rho,
			 &gas[i].Ne,
			 &gas[i].Neutral1H,
			 &gas[i].eps,
			 &gas[i].sfr);
		  gas[i].pos[0]=fac_long*gas[i].pos[0];
		  gas[i].pos[1]=fac_long*gas[i].pos[1];
		  gas[i].pos[2]=fac_long*gas[i].pos[2];
		  gas[i].vel[0]=fac_vel*gas[i].vel[0];
		  gas[i].vel[1]=fac_vel*gas[i].vel[1];
		  gas[i].vel[2]=fac_vel*gas[i].vel[2];
		  gas[i].masa=fac_masa*gas[i].masa;
		}
	    }
	}
      if(k>0)
	{
	  if(header.npartTotal[k]>0)
	    {
	      part[k-1]=malloc(header.npartTotal[k]*sizeof(struct nocolisional));
	      for(i=0;i<header.npartTotal[k];i++)
		{
		  fscanf(infile,"%d %f %f %f %f %f %f %f",
			 &part[k-1][i].id,
			 &part[k-1][i].pos[0],&part[k-1][i].pos[1],&part[k-1][i].pos[2],
			 &part[k-1][i].vel[0],&part[k-1][i].vel[1],&part[k-1][i].vel[2],
			 &part[k-1][i].masa);
		  part[k-1][i].pos[0]=fac_long*part[k-1][i].pos[0];
		  part[k-1][i].pos[1]=fac_long*part[k-1][i].pos[1];
		  part[k-1][i].pos[2]=fac_long*part[k-1][i].pos[2];
		  part[k-1][i].vel[0]=fac_vel*part[k-1][i].vel[0];
		  part[k-1][i].vel[1]=fac_vel*part[k-1][i].vel[1];
		  part[k-1][i].vel[2]=fac_vel*part[k-1][i].vel[2];
		  part[k-1][i].masa=fac_masa*part[k-1][i].masa;
		}
	    }
	}
      fclose(infile);
    }

  dummi=Npart*3*sizeof(float);
  printf("%d %d\n",dummi,Npart);
  //==========================================================
  //Write positios of the particles 
  //==========================================================
  fwrite(&dummi,sizeof(int),1,outfile);
  for(k=0;k<6;k++)
    {
      if(k==0)
	{
	  if(header.npartTotal[0]>0)
	    {
	      for(j=0;j<header.npartTotal[0];j++)
		{
		       
		      fwrite(&gas[j].pos,sizeof(float),3,outfile);		    
		}
	    }
	}
      if(k>0)
	{
	  if(header.npartTotal[k]>0)
	    {
	      for(j=0;j<header.npartTotal[k];j++)
		{
		  
		      fwrite(&part[k-1][j].pos,sizeof(float),3,outfile);
		    
		}
	    }
	}
    }
  fwrite(&dummi,sizeof(int),1,outfile);
  //==========================================================
  dummi=Npart*3*sizeof(float);
  printf("%d %d\n",dummi,Npart);
  //==========================================================
  //Write velocities of the particles 
  //==========================================================
  fwrite(&dummi,sizeof(int),1,outfile);
  for(k=0;k<6;k++)
    {
      if(k==0)
	{
	  if(header.npartTotal[0]>0)
	    {
	      for(j=0;j<header.npartTotal[0];j++)
		{
		  
		  
		      fwrite(&gas[j].vel,sizeof(float),3,outfile);
		    
		}
	    }
	}
      if(k>0)
	{
	  if(header.npartTotal[k]>0)
	    {
	      for(j=0;j<header.npartTotal[k];j++)
		{
		      fwrite(&part[k-1][j].vel,sizeof(float),3,outfile);
		}
	    }
	}
    }
  fwrite(&dummi,sizeof(int),1,outfile);
  //==========================================================

  dummi=Npart*sizeof(int);
  printf("%d %d\n",dummi,Npart);
  //==========================================================
  //Write id's of the particles 
  //==========================================================
  fwrite(&dummi,sizeof(int),1,outfile);
  for(k=0;k<6;k++)
    {
      if(k==0)
	{
	  if(header.npartTotal[0]>0)
	    {
	      for(j=0;j<header.npartTotal[0];j++)
		{
		  fwrite(&gas[j].id,sizeof(int),1,outfile);
		}
	    }
	}
      if(k>0)
	{
	  if(header.npartTotal[k]>0)
	    {
	      for(j=0;j<header.npartTotal[k];j++)
		{
		  fwrite(&part[k-1][j].id,sizeof(int),1,outfile);
		}
	    }
	}
    }
  fwrite(&dummi,sizeof(int),1,outfile);
  //==========================================================
  
  //==========================================================
  //Write masses of particles if necessary
  //==========================================================

  dummi=0;
  
  for(k=0;k<6;k++)
    {
      if((header.mass[k]==0.0) && (header.npartTotal[k]>0))
	{
	  dummi=dummi+header.npartTotal[k]*sizeof(float);
	}
    }

  printf("%d\n",dummi);
  
  if(dummi>0)
    {
      fwrite(&dummi,sizeof(int),1,outfile);
    }

  for(k=0;k<6;k++)
    {
      if(k==0)
	{
	  if((header.mass[0]==0.0) && (header.npartTotal[0]>0))
	    {
	      for(j=0;j<header.npartTotal[0];j++)
		{
		  fwrite(&gas[j].masa,sizeof(float),1,outfile);
		}
	    }
	  if((header.mass[0]!=0.0) && (header.npartTotal[0]>0))
	    {
	      printf("Las particulas de tipo %d tienen masa igual a %f\n",k,header.mass[0]);
	    }
	}
      if(k>0)
	{
	  
	  if((header.mass[k]==0.0) && (header.npartTotal[k]>0))
	    {
	      for(j=0;j<header.npartTotal[k];j++)
		{
		  fwrite(&part[k-1][j].masa,sizeof(float),1,outfile);
		}
	    }
	  if((header.mass[k]!=0.0) && (header.npartTotal[k]>0))
	    {
	      printf("Las particulas de tipo %d tienen masa igual a %f\n",k,header.mass[k]);
	    }
	}
      
    }
  if(dummi>0)
    {
      fwrite(&dummi,sizeof(int),1,outfile);
    }
  //==========================================================

  dummi=header.npartTotal[0]*sizeof(float);
  printf("%d\n",dummi);
  //Write properties of the particles of gas
  //==========================================================

  if(header.npartTotal[0]>0)
    {
      fwrite(&dummi,sizeof(int),1,outfile);
      for(j=0;j<header.npartTotal[0];j++)
	{
	  fwrite(&gas[j].U,sizeof(int),1,outfile);
	}
      fwrite(&dummi,sizeof(int),1,outfile);

      fwrite(&dummi,sizeof(int),1,outfile);
      for(j=0;j<header.npartTotal[0];j++)
	{
	  fwrite(&gas[j].rho,sizeof(int),1,outfile);
	}
      fwrite(&dummi,sizeof(int),1,outfile);

      fwrite(&dummi,sizeof(int),1,outfile);
      for(j=0;j<header.npartTotal[0];j++)
	{
	  fwrite(&gas[j].Ne,sizeof(int),1,outfile);
	}
      fwrite(&dummi,sizeof(int),1,outfile);

      fwrite(&dummi,sizeof(int),1,outfile);
      for(j=0;j<header.npartTotal[0];j++)
	{
	  fwrite(&gas[j].Neutral1H,sizeof(int),1,outfile);
	}
      fwrite(&dummi,sizeof(int),1,outfile);

      fwrite(&dummi,sizeof(int),1,outfile);
      for(j=0;j<header.npartTotal[0];j++)
	{
	  fwrite(&gas[j].eps,sizeof(int),1,outfile);
	}
      fwrite(&dummi,sizeof(int),1,outfile);

      fwrite(&dummi,sizeof(int),1,outfile);
      for(j=0;j<header.npartTotal[0];j++)
	{
	  fwrite(&gas[j].sfr,sizeof(int),1,outfile);
	}
      fwrite(&dummi,sizeof(int),1,outfile);
    }
  
  //==========================================================
  printf("%d\n",dummi);
  fclose(outfile);
}
