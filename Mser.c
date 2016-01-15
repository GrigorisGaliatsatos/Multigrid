#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include "mpi.h"

void Multigrid(double *L,double *U,double *eF,double *eB,double *b,double dx,double c,double d,double dt,double t,int M,int M_F,int N,int k,int n,FILE *fp3D);
void Gauss(double *X,double c,double d,double t,int M,int k,double *L,double *U,double *b);
void GSeidel(double *X,double *b,double c,double d,int M,int k,int MAX);
void Starting(double *L,double *U,double *X,double *b,double dx,double c,double d,double t,double dt,int M,int M_F,int N,int k,int n);
double normINF(double *X,int M,double T,double dx);
double norm2(double *X,int M,double T,double dx);
double f(double x,double t);
double u0(double t);
double u1(double t);
double V(double x);
double UAK(double x,double t);
void print_vector(double *X,int M);

int main(int argc,char **argv)
	{
	int i,M,N,k,n,tmp,rank,numtasks;
	double T,dx,dt,t,c,d,Time;
	double *X,*eF,*eB,*b,*L,*U;
	FILE *fp2D,*fp3D,*fpRES;
	MPI_Status stat;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);



//	fp3D=fopen("IIID_Schem.v","w");
//	fp2D=fopen("IID_Schem.v","w");
//	fpRES=fopen("Results.v","a");

	if(argc!=4)
		{
		printf("M=Ari8mos ypodiasthmatwn tou x.\n");
		printf("N=Ari8mos ypodiasthmatwn tou t.\n");
		printf("T=megistos xronos t, gia thn lysh.\n");
		printf("DWSTE:a.out [M] [N] [T]\n");
		MPI_Finalize();
		return 0;
		}
	else
		{
//		printf("DWSTE:[M] [N] [T]\n");
		M=atoi(*(argv+1));
		N=atoi(*(argv+2));
		T=atof(*(argv+3));
		}
	if((M<=0)||(N<=0)||(T<=0))
		{
		printf("La8ws eisodos dedomenwn !!\n");
		MPI_Finalize();
		return 0;
		}

	tmp=M;
	for(;;)
		{
		if(tmp==1)	break;
		if((int)(tmp%2)!=0)
			{
			printf("To M=Ari8mos ypodiasthmatwn tou x, prepei na einai pollaplasio tou 2.\n");
			MPI_Finalize();
			return 0;
			}
		tmp=tmp/2.0;
		}
		

//	Desmeysh mnumhs	:
	X=(double *)calloc(M+1,sizeof(double));
	eF=(double *)calloc(M+1,sizeof(double));
	eB=(double *)calloc(M+1,sizeof(double));
	b=(double *)calloc(M+1,sizeof(double));
	L=(double *)calloc(M-1,sizeof(double));
	U=(double *)calloc(M,sizeof(double));

//	Enarksh metrhshs xronou ekteleshs praksewn	:
	Time=((double)(clock())/CLOCKS_PER_SEC);

//	Arxikopoihsh timwn	:
	dx=1.0/(double)(M);
	dt=T/(double)(N);
	k=1;
	c=((dx*dx)+2.0*dt)/dt;
	d=-1.0;
	for (i=0;i<=M;i++)
		{
		*(X+i)=V((double)(i)*dx);
		*(eF+i)=0.0;
		*(eB+i)=0.0;
		}
//	Lush diaforikhs	:

//	for(i=0;i<M+1;i++)	fprintf(fp3D,"%.6lf %.6lf %.6lf %.6lf\n",((double)(i)*dx),0.0,*(X+i),UAK((double)(i)*dx,0.0));

	if(M<=64)
		{
		for(n=1;n<=N;n++)
			{
			t=(double)(n)*dt;
			for(i=1;i<M;i++)	*(b+i)=(dx*dx)*( (*(X+i)/dt) + f((double)(i)*dx,t) );
			*(b+1)=*(b+1)+u0(t);
			*(b+(M-1))=*(b+(M-1))+u1(t);
			Gauss(X,c,d,t,M,1,L,U,b);
//			for(i=0;i<=M;i++)	fprintf(fp3D,"%.6lf %.6lf %.6lf %.6lf\n",((double)(i)*dx),t,*(X+i),UAK((double)(i)*dx,t));
			for(i=0;i<=M;i++)	printf("%.6lf %.6lf %.6lf %.6lf\n",((double)(i)*dx),t,*(X+i),UAK((double)(i)*dx,t));
			}
		}
	else
		{
		for(n=1;n<=N;n++)
			{
			t=(double)(n)*dt;
			for(i=1;i<=M-1;i++)	*(b+i)=(dx*dx)*((*(X+i)/dt)+f((double)(i)*dx,t));
			*(b+1)=*(b+1)+u0(t);
			*(b+(M-1))=*(b+(M-1))+u1(t);
			Starting(L,U,X,b,dx,c,d,t,dt,M,M,N,k,n);

			GSeidel(X,b,c,d,M,k,3);

			*(b+1)=*(b+1)-c*(*(X+1))-d*(*(X+2));
			for(i=2;i<=M-2;i++)	*(b+i)=*(b+i)-d*(*(X+i-1))-c*(*(X+i))-d*(*(X+i+1));
			*(b+M-1)=*(b+M-1)-d*(*(X+M-2))-c*(*(X+M-1));

//			printf("\n To dianusma r=b-Au einai to ekshs:\n");
//			print_vector(b,M);

			Multigrid(L,U,eF,eB,b,dx,c,d,dt,t,M,M,N,k,n,fp3D);
			for (i=0;i<=M;i++)	*(eF+i)=0.0;

//			printf("\n To dianusma eB sto Xroniko bhma n=%d, einai to ekshs:\n",n);
//			print_vector(eB,M);

			for(i=1;i<=M-1;i++)	*(X+i)=*(X+i)+(*(eB+i));
//			for(i=0;i<=M;i++)	fprintf(fp3D,"%.6lf %.6lf %.6lf %.6lf\n",((double)(i)*dx),t,*(X+i),UAK((double)(i)*dx,t));
			}
		}
	Time=((double)(clock())/CLOCKS_PER_SEC)-Time;

//	for(i=0;i<M+1;i++)	fprintf(fp2D,"%.6lf %.6lf %.6lf %.6lf\n",((double)(i)*dx),T,*(X+i),UAK((double)(i)*dx,T));
//	fprintf(fpRES,"M=%.6d N=%.6d T=%.6lf Sfalma=%.6lf Time=%.6lf\n",M,N,T,norm2(X,M,T,dx),Time);
	printf("\nM=%.6d N=%.6d T=%.6lf Norm2=%.6lf NormINF=%.6lf L^2=%.6lf Time=%.6lf\n",M,N,T,norm2(X,M,T,dx),normINF(X,M,T,dx),norm2(X,M,T,dx)*sqrt(dx),Time);

//	fclose(fpRES);
//	fclose(fp2D);
//	fclose(fp3D);
	free(L);
	free(U);
	free(b);
	free(X);
	free(eF);
	free(eB);

	MPI_Finalize();

	return 0;
	}

void Multigrid(double *L,double *U,double *eF,double *eB,double *b,double dx,double c,double d,double dt,double t,int M,int M_F,int N,int k,int n,FILE *fp3D)
	{
	int i;
	double *v,*r;

	if(k!=1)
		{
		v=(double *)calloc(M,sizeof(double));
		r=(double *)calloc(M,sizeof(double));
		for(i=1;i<=M-1;i++)	*(v+i)=(*(b+i*2-1)+2.0*(*(b+2*i))+*(b+2*i+1))/4.0;
		}

	if(((dx>(3.0)*(1.0/(double)(M_F)))&&(500.0>((1.0/dx)+1.0)))||(k>=8))
		{
		for (i=1;i<=M_F-1;i++)	*(eB+i)=0.0;
		Gauss(eB,c,d,t,M,k,L,U,v);
		*eB=0.0;
		*(eB+M_F)=0.0;	
		}
	else
		{
		if(k!=1)
			{
			GSeidel(eF,v,c,d,M,k,3);

//			printf("\n To dianusma e[%d] einai to ekshs:\n",k);
//			print_vector(eF,M);
			for(i=0;i<=M-1;i++)	*(r+i)=*(v+i);

//			printf("\n To dianusma r[%d]=b-Au einai to ekshs:\n",k);
//			print_vector(r,M);

			*(v+1)=*(v+1)-c*(*(eF+k))-d*(*(eF+2*k));
			for(i=2;i<=M-2;i++)	*(v+i)=*(v+i)-d*(*(eF+(i-1)*k))-c*(*(eF+i*k))-d*(*(eF+(i+1)*k));
			*(v+M-1)=*(v+M-1)-d*(*(eF+(M-2)*k))-c*(*(eF+(M-1))*k);
			}

//		printf("\n To dianusma r sto bhma k=%d, einai to ekshs:\n",k);
//		if(k==1)	print_vector(b,M);
//		else	print_vector(v,M);

		if(k==1)	Multigrid(L,U,eF,eB,b,2.0*dx,(6.0*c+8.0*d)/8.0,(c+4.0*d)/8.0,dt,t,(int)(M/2),M_F,N,(int)(2*k),n,fp3D);
		else	Multigrid(L,U,eF,eB,v,2.0*dx,(6.0*c+8.0*d)/8.0,(c+4.0*d)/8.0,dt,t,(int)(M/2),M_F,N,(int)(2*k),n,fp3D);

//		Grammikh parembolh	:
		for(i=1;i<=M-1;i+=2)	*(eB+i*k)=(*(eB+(i-1)*k)+*(eB+(i+1)*k))/2.0;

//		Klhsh ths Gauss-Seidel gia thn anodo tou V-circle	:
		if(k==1)	GSeidel(eB,b,c,d,M,k,3);
		if(k!=1)	GSeidel(eB,r,c,d,M,k,3);
		}
	if(k!=1)
		{
		free(v);
		free(r);
		}
	}

void Starting(double *L,double *U,double *X,double *b,double dx,double c,double d,double t,double dt,int M,int M_F,int N,int k,int n)
	{
	int i;
	double *v;

	if(k!=1)
		{
		v=(double *)calloc(M,sizeof(double));
		for(i=1;i<=M-1;i++) *(v+i)=(*(b+i*2-1)+2.0*(*(b+2*i))+*(b+2*i+1))/4.0;
		}

	if(((dx>(3.0)*(1.0/(double)(M_F)))&&(500.0>((1.0/dx)+1.0)))||(k>=8))
		{
		Gauss(X,c,d,t,M,k,L,U,v);
		}
	else
		{
		if(k==1)	Starting(L,U,X,b,2.0*dx,(6.0*c+8.0*d)/8.0,(c+4.0*d)/8.0,t,dt,(int)(M/2),M_F,N,(int)(2*k),n);
		else	Starting(L,U,X,v,2.0*dx,(6.0*c+8.0*d)/8.0,(c+4.0*d)/8.0,t,dt,(int)(M/2),M_F,N,(int)(2*k),n);

//		Grammikh parembolh	:
		for(i=1;i<=M-1;i+=2)	*(X+i*k)=(*(X+(i-1)*k)+*(X+(i+1)*k))/2.0;
		}
	if(k!=1) free(v);
	}

void Gauss(double *X,double c,double d,double t,int M,int k,double *L,double *U,double *b)
	{
	int i;
	double *y;

	y=(double *)calloc(M-1,sizeof(double));

//	L-U_decomposition	:
	*U=c;
	*L=d/(*U);
	for(i=0;i<M-1;i++)
		{
		if((i>0)&&(i<M-2))	*(L+i)=d/(*(U+i));
		*(U+i+1)=c-d*(*(L+i));
		}

//	Solve_Lower_Triangular	:
	*(y+0)=*(b+1);
	for(i=1;i<M-1;i++)	*(y+i)=(*(b+(i+1))-(*(L+i-1))*(*(y+i-1)));

//	solve_upper_triangular	:
	*(X+(M-1)*k)=(*(y+M-2))/(*(U+M-2));
	for(i=M-3;i>=0;i--)	*(X+(i+1)*k)=(*(y+i)-d*(*(X+(i+2)*k)))/(*(U+i));
	*X=u0(t);
	*(X+M*k)=u1(t);

//	printf("\n Ta metaballomena stoixeia tou L einai ta ekshs:\n");
//	print_vector(L,M-3);
//	printf("\n Ta metaballomena stoixeia tou U einai ta ekshs:\n");
//	print_vector(U,M-2);
//	printf("\n To dianusma b einai to ekshs:\n");
//	print_vector(b,M);
//	printf("\n To dianusma y (Pou einai h lush tou susthmatos L*y=b) einai to ekshs:\n");
//	print_vector(y,M-2);
//	printf("\n To dianusma X (Pou einai h lush tou susthmatos A*X=b <=> L*U*X=b) einai to ekshs:\n");
//	print_vector(X,M);

	free(y);
	}

void GSeidel(double *x,double *b,double c,double d,int M,int k,int MAX)
	{
	double *y;
	int i,j;

	y=(double *)malloc(M*sizeof(double ));

	for(i=1;i<=M-1;i++)	*(y+i)=*(x+i*k);
	for(j=0;j<MAX;j++)
		{
		*(x+k)=(*(b+1)-d*(*(y+2)))/c;
		for(i=1;i<=M-2;i++)	*(x+(i+1)*k)=(-1.0*d*(*(x+i*k))-d*(*(y+(i+2)))+(*(b+i+1)))/c;
		*(x+(M-1)*k)=((*(b+(M-1)))-d*(*(x+(M-2)*k)))/c;
		for(i=1;i<=M-1;i++)	*(y+i)=*(x+i*k);
		}
	free(y);
	}

double normINF(double *X,int M,double T,double dx)
	{
	int i;
	double tmp,nrINF=0.0;

	for(i=0;i<=M;i++)
		{
		tmp=*(X+i)-UAK((double)(i)*dx,T);
		if(tmp<0)	tmp=-tmp;
		if(nrINF<tmp)	nrINF=tmp;
		}

	return nrINF;
	}

double norm2(double *X,int M,double T,double dx)
	{
	int i;
	double tmp=0.0,nr2=0.0;

	for(i=0;i<=M;i++)
		{
		tmp=tmp+(*(X+i)-UAK((double)(i)*dx,T))*(*(X+i)-UAK((double)(i)*dx,T));
		
		}
        nr2=sqrt(tmp);
	return nr2;
	}

void print_vector(double *X,int M)
	{
	int i;
	for(i=0;i<=M;i++)
		{
		if(*(X+i)<0.0)	printf("%.12lf\n",*(X+i));
		else	printf(" %.12lf\n",*(X+i));
		}
	}

double f(double x,double t)
	{
	double F;
//	F=(3.0*t*t)-(6.0*x);
//	F=cos(t)+sin(x);
	F=(16.0*cos(16.0*t))+(16.0*16.0*(sin(16.0*x)));
//	F=cos(t)+(1.0-2.0*x)*(1.0-2.0*x)*sin(x*(1.0-x))+2.0*cos(x*(1.0-x));
	return F;
	}

double u0(double t)
	{
	double U0;
//	U0=t*t*t;
//	U0=sin(t);
	U0=sin(t*16.0);
//	U0=sin(t);
	return U0;
	}

double u1(double t)
	{
	double U1;
//	U1=(t*t*t)+1.0;
//	U1=sin(1.0)+sin(t);
	U1=sin(16.0)+sin(t*16.0);
//	U1=sin(t);
	return U1;
	}

double V(double x)
	{
	double v;
//	v=x*x*x;
//	v=sin(x);
	v=sin(16.0*x);
//	v=sin(x);
	return v;
	}

double UAK(double x,double t)
	{
	double AK;
//	AK=(x*x*x)+(t*t*t);
//	AK=sin(x)+sin(t);
	AK=sin(16.0*x)+sin(16.0*t);
//	AK=sin(x*(1.0-x))+sin(t);
	return AK;
	}

