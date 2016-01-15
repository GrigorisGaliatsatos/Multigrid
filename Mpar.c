#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include "mpi.h"

void Multigrid(double *L,double *U,double *eF,double *eB,double *b,double *btmp,double dx,double c,double d,double dt,double t,int M,int M_F,int N,int k,int n,int rank,int numtasks,MPI_Status stat,FILE *fp3D);
void Gauss(double *X,double c,double d,double t,int M,int k,double *L,double *U,double *b);
void GSeidel(double *X,double *b,double *btmp,double c,double d,int M,int k,int MAX,int rank,int numtasks,MPI_Status stat);
void Starting(double *L,double *U,double *X,double *b,double *btmp,double dx,double c,double d,double t,double dt,int M,int M_F,int N,int k,int n,int rank,int numtasks,MPI_Status stat);
double normINF(double *X,int M,double T,double dx,int rank);
double Spow2(double *X,int M,double T,double dx,int rank);
double f(double x,double t);
double u0(double t);
double u1(double t);
double V(double x);
double UAK(double x,double t);
void print_vector(double *X,int M);

int main(int argc,char **argv)
	{
	int i,M,Mp,N,k,n,tmp,rank,numtasks;
	double dtmp,T,dx,dt,t,c,d,norm2,nrmINF,Time;
	double *X,*eF,*eB,*b,*btmp,*L,*U;
	FILE *fp2D,*fp3D,*fpRES;
	MPI_Status stat;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if(rank==0)
		{
//		fp3D=fopen("IIID_Schem.v","w");
//		fp2D=fopen("IID_Schem.v","w");
//		fpRES=fopen("Results.v","a");
		if(argc!=4)
			{
			printf("M=Ari8mos ypodiasthmatwn tou x.\n");
			printf("N=Ari8mos ypodiasthmatwn tou t.\n");
			printf("T=megistos xronos t, gia thn lysh.\n");
			printf("DWSTE:a.out [M] [N] [T]\n");
			return 0;
			}
		else
			{
			M=atoi(*(argv+1));
			N=atoi(*(argv+2));
			T=atof(*(argv+3));
			}
		if((M<=64)||(N<=0)||(T<=0))
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
			         printf("To M=Ari8mos ypodiasthmatwn tou x, prepei na einai dynamh tou 2.\n");
				 MPI_Finalize();
        			 return 0;
				}
			tmp=tmp/2.0;
			}
		
		}

	MPI_Bcast(&M,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&T,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

//	Enarksh metrhshs xronou ekteleshs praksewn	:
	Time=((double)(clock())/CLOCKS_PER_SEC);

//	Arxikopoihsh timwn	:
	Mp=M/numtasks;

//	Desmeysh mnumhs	:
	X=(double *)calloc(Mp+2,sizeof(double));

   
	eF=(double *)calloc(Mp+2,sizeof(double));
	eB=(double *)calloc(Mp+1,sizeof(double));
	b=(double *)calloc(Mp+1,sizeof(double));
	if(rank==0)	btmp=(double *)calloc(numtasks,sizeof(double));
	if(rank==numtasks-1)
		{
		L=(double *)calloc(M,sizeof(double));
		U=(double *)calloc(M+1,sizeof(double));
	   	}
//	Arxikopoihsh timwn	:
	dx=1.0/(double)(M);
	dt=T/(double)(N);
	k=1;
	c=((dx*dx)+2.0*dt)/dt;
	d=-1.0;
	if(rank==0) *X=V(0.0);
	for (i=1;i<=Mp;i++)
		{
		*(X+i)=V((double)(i+rank*Mp)*dx);
		*(eF+i)=0.0;
		*(eB+i)=0.0;
		}

//	Lush diaforikhs	:

//	for(i=0;i<M+1;i++)	fprintf(fp3D,"%.6lf %.6lf %.6lf %.6lf\n",((double)(i)*dx),0.0,*(X+i),UAK((double)(i)*dx,0.0));

	for(n=1;n<=N;n++)
		{
//		printf("Briskomai sto xroniko bhma %d, ek twn %d\n",n,N);
		t=(double)(n)*dt;

		if(rank!=numtasks-1) for(i=1;i<=Mp;i++)	*(b+i)=(dx*dx)*((*(X+i)/dt)+f((double)(i+rank*Mp)*dx,t));

		if(rank==numtasks-1) for(i=1;i<=Mp-1;i++)	*(b+i)=(dx*dx)*((*(X+i)/dt)+f((double)(i+rank*Mp)*dx,t));

		if(rank==0) *(b+1)=*(b+1)+u0(t);

		if(rank==numtasks-1) *(b+(Mp-1))=*(b+(Mp-1))+u1(t);
		

      		if(rank!=0) MPI_Send(b+1,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
		else for(i=1;i<=numtasks-1;i++) MPI_Recv(btmp+i-1,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&stat);

		if(rank==0)
			{
			for(i=1;i<=numtasks-1;i++) MPI_Send(btmp+i,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD);
			*b=*btmp;
			}
		else MPI_Recv(b,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,&stat);

		Starting(L,U,X,b,btmp,dx,c,d,t,dt,Mp+1,M,N,k,n,rank,numtasks,stat);

		GSeidel(X,b,btmp,c,d,Mp+1,k,3,rank,numtasks,stat);

		if(rank!=0) MPI_Send(X+1,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
		else for(i=1;i<=numtasks-1;i++) MPI_Recv(btmp+i-1,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&stat);

		if(rank==0)
			{
			for(i=1;i<=numtasks-1;i++) MPI_Send(btmp+i,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD);
			*(X+Mp+1)=*btmp;
			}
		else MPI_Recv(X+Mp+1,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,&stat);

		if(rank==0)	*(b+1)=*(b+1)-c*(*(X+1))-d*(*(X+2));
   		else	*(b+1)=*(b+1)-d*(*X)-c*(*(X+1))-d*(*(X+2));
		for(i=2;i<=Mp-2;i++)	*(b+i)=*(b+i)-d*(*(X+i-1))-c*(*(X+i))-d*(*(X+i+1));
		if(rank!=numtasks-1)
			{
			*(b+Mp-1)=*(b+Mp-1)-d*(*(X+Mp-2))-c*(*(X+Mp-1))-d*(*(X+Mp));
			*(b+Mp)=*(b+Mp)-d*(*(X+Mp-1))-c*(*(X+Mp))-d*(*(X+Mp+1));
			}
		else	*(b+Mp-1)=*(b+Mp-1)-d*(*(X+Mp-2))-c*(*(X+Mp-1));

//		printf("\n To dianusma r=b-Au einai to ekshs:\n");
//		print_vector(b,M);

		Multigrid(L,U,eF,eB,b,btmp,dx,c,d,dt,t,Mp+1,M,N,k,n,rank,numtasks,stat,fp3D);
		for (i=0;i<=Mp;i++)	*(eF+i)=0.0;

//		printf("\n To dianusma eB sto Xroniko bhma n=%d, einai to ekshs:\n",n);
//		print_vector(eB,M);

		for(i=1;i<=Mp;i++)	*(X+i)=*(X+i)+(*(eB+i));
//		for(i=0;i<=M;i++)	fprintf(fp3D,"%.6lf %.6lf %.6lf %.6lf\n",((double)(i)*dx),t,*(X+i),UAK((double)(i)*dx,t));
		}

	Time=((double)(clock())/CLOCKS_PER_SEC)-Time;
	dtmp=Spow2(X,Mp,T,dx,rank);
	MPI_Reduce(&dtmp,&norm2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	dtmp=normINF(X,Mp,T,dx,rank);
	MPI_Reduce(&dtmp,&nrmINF,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
//	for(i=0;i<M+1;i++)	fprintf(fp2D,"%.6lf %.6lf %.6lf %.6lf\n",((double)(i)*dx),T,*(X+i),UAK((double)(i)*dx,T));
// 	fprintf(fpRES,"M=%.6d N=%.6d T=%.6lf Sfalma=%.6lf Time=%.6lf\n",M,N,T,norm2(X,M,T,dx),Time);
	if(rank==0)	printf("\nProc=%d M=%.6d N=%.6d T=%.6lf Norm2=%.6lf NormINF=%.6lf L^2=%.6lf Time=%.6lf\n",numtasks,M,N,T,sqrt(norm2),nrmINF,sqrt(norm2*dx),Time);

//	fclose(fpRES);
//	fclose(fp2D);
//	fclose(fp3D);
	if(rank==numtasks-1)
		{
		free(L);
		free(U);
		}
	free(b);
        if(rank==0)	free(btmp);
	free(X);
	free(eF);
	free(eB);

	MPI_Finalize();


	return 0;
	}

void Multigrid(double *L,double *U,double *eF,double *eB,double *b,double *btmp,double dx,double c,double d,double dt,double t,int M,int M_F,int N,int k,int n,int rank,int numtasks,MPI_Status stat,FILE *fp3D)
	{
	int i;
	double *v,*r,*vB,*eB_B,tmpeF;

	if(k!=1)
		{
		v=(double *)calloc(M,sizeof(double));
		r=(double *)calloc(M,sizeof(double));

		if(rank!=0) MPI_Send(b+1,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
		else for(i=1;i<=numtasks-1;i++) MPI_Recv(btmp+i-1,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&stat);

		if(rank==0)
			{
			for(i=1;i<=numtasks-1;i++) MPI_Send(btmp+i,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD);
			*b=*btmp;
			}
		else MPI_Recv(b,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,&stat);

		for(i=1;i<=M-2;i++)	*(v+i)=(*(b+i*2-1)+2.0*(*(b+2*i))+*(b+2*i+1))/4.0;
		if(rank!=numtasks-1) *(v+M-1)=(*(b+2*(M-1)-1)+2.0*(*(b+2*(M-1)))+*b)/4.0;
		}

	if(((dx>(3.0)*(1.0/(double)(M_F)))&&(500.0>((1.0/dx)+1.0)))||(k>=8))
		{
	for(i=0;i<=(int)(M_F/numtasks);i++) *(eB+i)=0.0;

	if(rank!=numtasks-1) MPI_Send(v+1,M-1,MPI_DOUBLE,numtasks-1,rank,MPI_COMM_WORLD);
	else
		{
		eB_B=(double *)calloc((M-1)*numtasks+1,sizeof(double));
		vB=(double *)calloc((M-1)*numtasks,sizeof(double));

		for (i=1;i<=(M-1)*numtasks-1;i++)	*(eB_B+i)=0.0;

		for(i=0;i<=numtasks-2;i++) MPI_Recv(vB+1+i*(M-1),M-1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&stat);
		for(i=1;i<=M-2;i++)	*(vB+(numtasks-1)*(M-1)+i)=*(v+i);
		Gauss(eB_B,c,d,t,(M-1)*numtasks,1,L,U,vB);
		*eB_B=0.0;
		*(eB_B+(M-1)*numtasks)=0.0;
		for(i=0;i<=numtasks-2;i++)
		MPI_Send(eB_B+i*(M-1),M,MPI_DOUBLE,i,i,MPI_COMM_WORLD);

		free(vB);
		}
		if(rank!=numtasks-1)
			{
			eB_B=(double *)calloc(M,sizeof(double));
			MPI_Recv(eB_B,M,MPI_DOUBLE,numtasks-1,rank,MPI_COMM_WORLD,&stat);
			for(i=0;i<=M-1;i++) *(eB+i*k)=*(eB_B+i);
			free(eB_B);
			}
		else
			{
			for(i=0;i<=M-1;i++) *(eB+i*k)=*(eB_B+(M-1)*(numtasks-1)+i);
			free(eB_B);
			}
		}
	else
		{
		if(k!=1)
			{
			GSeidel(eF,v,btmp,c,d,M,k,3,rank,numtasks,stat);
//			printf("\n To dianusma e[%d] einai to ekshs:\n",k);
//			print_vector(eF,M);
			for(i=0;i<=M-1;i++)	*(r+i)=*(v+i);

//			printf("\n To dianusma r[%d]=b-Au einai to ekshs:\n",k);
//			print_vector(r,M);


		if(rank!=0) MPI_Send(eF+k,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
		else for(i=1;i<=numtasks-1;i++) MPI_Recv(btmp+i-1,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&stat);

		if(rank==0)
			{
			for(i=1;i<=numtasks-1;i++) MPI_Send(btmp+i,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD);
			tmpeF=*btmp;
			}
		else MPI_Recv(&tmpeF,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,&stat);

		if(rank==0)	*(v+1)=*(v+1)-c*(*(eF+k))-d*(*(eF+2*k));
		else	*(v+1)=*(v+1)-d*(*eF)-c*(*(eF+k))-d*(*(eF+2*k));

		for(i=2;i<=M-3;i++)	*(v+i)=*(v+i)-d*(*(eF+(i-1)*k))-c*(*(eF+i*k))-d*(*(eF+(i+1)*k));

		if(rank!=numtasks-1)
			{
			*(v+M-2)=*(v+M-2)-d*(*(eF+(M-3)*k))-c*(*(eF+(M-2)*k))-d*(*(eF+(M-1)*k));
			*(v+M-1)=*(v+M-1)-d*(*(eF+(M-2)*k))-c*(*(eF+(M-1)*k))-d*(tmpeF);
			}
		else	*(v+M-2)=*(v+M-2)-d*(*(eF+(M-3)*k))-c*(*(eF+(M-2)*k));
		}

//		printf("\n To dianusma r sto bhma k=%d, einai to ekshs:\n",k);
//		if(k==1)	print_vector(b,M);
//		else	print_vector(v,M);

		if(k==1)	Multigrid(L,U,eF,eB,b,btmp,2.0*dx,(6.0*c+8.0*d)/8.0,(c+4.0*d)/8.0,dt,t,(int)((M-1)/2)+1,M_F,N,(int)(2*k),n,rank,numtasks,stat,fp3D);
		else	Multigrid(L,U,eF,eB,v,btmp,2.0*dx,(6.0*c+8.0*d)/8.0,(c+4.0*d)/8.0,dt,t,(int)((M-1)/2)+1,M_F,N,(int)(2*k),n,rank,numtasks,stat,fp3D);

//		Grammikh parembolh	:
		for(i=1;i<=M-1;i+=2)	*(eB+i*k)=(*(eB+(i-1)*k)+(*(eB+(i+1)*k)))/2.0;

//		Klhsh ths Gauss-Seidel gia thn anodo tou V-circle	:
		if(k==1)	GSeidel(eB,b,btmp,c,d,M,k,3,rank,numtasks,stat);
		if(k!=1)	GSeidel(eB,r,btmp,c,d,M,k,3,rank,numtasks,stat);
		}
	if(k!=1)
   		{
		free(v);
		free(r);
		}
	}

void Starting(double *L,double *U,double *X,double *b,double *btmp,double dx,double c,double d,double t,double dt,int M,int M_F,int N,int k,int n,int rank,int numtasks,MPI_Status stat)
	{
	int i;
	double *v,*XB,*vB;

	if(k!=1)
		{
		v=(double *)calloc(M,sizeof(double));

		if(rank!=0) MPI_Send(b+1,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
		else for(i=1;i<=numtasks-1;i++) MPI_Recv(btmp+i-1,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&stat);

		if(rank==0)
			{
			for(i=1;i<=numtasks-1;i++) MPI_Send(btmp+i,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD);
			*b=*btmp;
			}
		else MPI_Recv(b,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,&stat);

		for(i=1;i<=M-2;i++)	*(v+i)=(*(b+i*2-1)+2.0*(*(b+2*i))+*(b+2*i+1))/4.0;
		if(rank!=numtasks-1) *(v+M-1)=(*(b+2*(M-1)-1)+2.0*(*(b+2*(M-1)))+*b)/4.0;
		}

	if(((dx>(3.0)*(1.0/(double)(M_F)))&&(500.0>((1.0/dx)+1.0)))||(k>=8))
		{

		if(rank!=numtasks-1) MPI_Send(v+1,M-1,MPI_DOUBLE,numtasks-1,rank,MPI_COMM_WORLD);
		else
			{
			XB=(double *)calloc((M-1)*numtasks+1,sizeof(double));
			vB=(double *)calloc((M-1)*numtasks,sizeof(double));

			for(i=0;i<=numtasks-2;i++) MPI_Recv(vB+1+i*(M-1),M-1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&stat);
			for(i=1;i<=M-2;i++)	*(vB+(numtasks-1)*(M-1)+i)=*(v+i);
			Gauss(XB,c,d,t,(M-1)*numtasks,1,L,U,vB);
			for(i=0;i<=numtasks-2;i++)
			MPI_Send(XB+i*(M-1),M,MPI_DOUBLE,i,i,MPI_COMM_WORLD);

			free(vB);
			}
		if(rank!=numtasks-1)
			{
			XB=(double *)calloc(M,sizeof(double));
			MPI_Recv(XB,M,MPI_DOUBLE,numtasks-1,rank,MPI_COMM_WORLD,&stat);
			for(i=0;i<=M-1;i++) *(X+i*k)=*(XB+i);
			free(XB);
			}
		else
			{
			for(i=0;i<=M-1;i++) *(X+i*k)=*(XB+(M-1)*(numtasks-1)+i);
			free(XB);
			}
		}
	else
		{
		if(k==1)	Starting(L,U,X,b,btmp,2.0*dx,(6.0*c+8.0*d)/8.0,(c+4.0*d)/8.0,t,dt,(int)((M-1)/2)+1,M_F,N,(int)(2*k),n,rank,numtasks,stat);
		else	Starting(L,U,X,v,btmp,2.0*dx,(6.0*c+8.0*d)/8.0,(c+4.0*d)/8.0,t,dt,(int)((M-1)/2)+1,M_F,N,(int)(2*k),n,rank,numtasks,stat);

//		Grammikh parembolh	:
		for(i=1;i<=M-1;i+=2)	*(X+i*k)=(*(X+(i-1)*k)+*(X+(i+1)*k))/2.0;
//		if(rank!=numtasks-1)	*(X+(M-1)*k)=(*(X+(M-2)*k)+*(X+M*k))/2.0;
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

void GSeidel(double *x,double *b,double *btmp,double c,double d,int M,int k,int MAX,int rank,int numtasks,MPI_Status stat)
	{
	double *y,xtmp;
	int i,j;

	y=(double *)malloc(M*sizeof(double ));

	for(i=2;i<=M-3;i+=2)	*(y+i)=*(x+i*k);
	if(rank!=numtasks-1)	*(y+M-1)=*(x+(M-1)*k);

	for(j=0;j<MAX;j++)
		{
		if((rank!=0)&&(rank!=numtasks-1))	MPI_Send(y+M-1,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
		if(rank==0)
			{
			for(i=1;i<=numtasks-2;i++)	MPI_Recv(btmp+i,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&stat);
			*btmp=*(y+M-1);
			}
		if(rank==0)	for(i=1;i<=numtasks-1;i++)	MPI_Send(btmp+i-1,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD);
		else	MPI_Recv(y,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,&stat);

		if(rank==0)	*(x+k)=(*(b+1)-d*(*(y+2)))/c;
		else	*(x+k)=(*(b+1)-d*(*(y+2))-d*(*y))/c;

		for(i=3;i<=M-4;i+=2)	*(x+i*k)=(-1.0*d*(*(y+i-1))-d*(*(y+i+1))+(*(b+i)))/c;

		if(rank==numtasks-1)	*(x+(M-2)*k)=((*(b+(M-2)))-d*(*(y+(M-3))))/c;
		else	*(x+(M-2)*k)=((*(b+M-2))-d*(*(y+M-3))-d*(*(y+M-1)))/c;
		
		if(rank!=0)	MPI_Send(x+k,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
		else	for(i=1;i<=numtasks-1;i++)	MPI_Recv(btmp+i-1,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&stat);

		if(rank==0)
			{
			xtmp=*btmp;
			for(i=1;i<=numtasks-1;i++)	MPI_Send(btmp+i,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD);
			}
		else	MPI_Recv(&xtmp,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,&stat);

		for(i=2;i<=M-3;i+=2)	*(x+i*k)=(-1.0*d*(*(x+(i-1)*k))-d*(*(x+(i+1)*k))+(*(b+i)))/c;
		if(rank!=numtasks-1) 	*(x+(M-1)*k)=(-1.0*d*(*(x+(M-2)*k))-d*xtmp+(*(b+M-1)))/c;

		for(i=2;i<=M-3;i+=2)	*(y+i)=*(x+i*k);
		if(rank!=numtasks-1)	*(y+M-1)=*(x+(M-1)*k);
		}
	free(y);
	}

double normINF(double *X,int M,double T,double dx,int rank)
	{
	int i;
	double tmp,nrINF=0.0;

	for(i=1;i<=M;i++)
		{
		tmp=*(X+i)-UAK((double)(i+rank*M)*dx,T);
		if(tmp<0)	tmp=-tmp;
		if(nrINF<tmp)	nrINF=tmp;
		}

	return nrINF;
	}

double Spow2(double *X,int M,double T,double dx,int rank)
	{
	int i;
	double tmp=0.0;

	for(i=1;i<=M;i++)	tmp=tmp+(*(X+i)-UAK((double)(i+rank*M)*dx,T))*(*(X+i)-UAK((double)(i+rank*M)*dx,T));
	return tmp;
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

