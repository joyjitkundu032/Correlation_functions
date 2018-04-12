#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "ran2.c"

#define eps 0.0000000001
#define D 6
#define D2 (3*3*3*3*3*3)
#define M 5
#define N 2000
#define Ravg 1.0
#define swapprob 0.2
#define pi (22.0/7.0)
#define tnum 100
#define dr 0.01
#define bin (1.0/dr)
#define Knum 8
#define qbin 4.2
#define SAMP 50000000
#define FRAC 0.20
#define bingof1 1000
#define buffer 1.05

int RMAX,RMAX_in,GAP2,nsteps,tstart,Kdim,Mcell,dim;
double boxsize,boxl,boxl2,rd,Rmin,Rdif,Rbound,vf,Kmin,cellsize,Rskin;
double pos[D][N],RD[N],rd_sum,count_in,numaa;
char readfile[200],outfile1[200],outfile2[200],outfile3[200];
double *count,*count_out,*num,*Sq_re,*Sqaa_re;
double *Sq_im,*Sqaa_im,*Sq,*Sqaa,*Sqbb;
double *g_of_1,*nwg_of_1;
long int seed=485620;

void read_input(int tt)
{
        int i,time,nci[D];
        double xi[D];
        time=tstart+tt*GAP2;
        FILE *fpr;
        sprintf(readfile,"./rho_%1.4lfSWP_%1.2lfeq/config%dD_N%dBS%1.4lfRmax%1.3lfRmin%1.3lf_vf%1.4lfRskin_%1.2lft_%d.dat",vf,swapprob,D,N,boxsize,rd,Rmin,vf,Rskin,time);
        fpr=fopen(readfile,"r");
        i=0;
        while(fscanf(fpr,"%lf%lf%lf%lf%lf%lf%d%d%d%d%d%d%lf%lf%lf%lf%lf%lf%lf",&pos[0][i],&pos[1][i],&pos[2][i],&pos[3][i],&pos[4][i],&pos[5][i],&nci[0],&nci[1],&nci[2],&nci[3],&nci[4],&nci[5],&xi[0],&xi[1],&xi[2],&xi[3],&xi[4],&xi[5],&RD[i])!=EOF)
                i++;
        fclose(fpr);
}

double anint(double x)
{
	double d;
	if(x>=0.50)
		d=1.00;
	else
	{
		if(x<=(-0.50))
			d=-1.00;
		else
			d=0.00;
	}
	return d;
}

void calculate_g_of_r()
{
	int k,i,j,rr,srr;
	double sep[D],sepsq,sepsqrt,rij,srij;
	double fac;	
	rd_sum=0.0;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			if(i !=j)
			{
				sepsq=0.0;
				for(k=0;k<D;k++)
				{
        	                        sep[k]=pos[k][i]-pos[k][j];
        	                        sep[k]=sep[k]-(boxl)*anint(sep[k]/(boxl));
        	                        sepsq=sepsq+sep[k]*sep[k];
				}
				sepsqrt=sqrt(sepsq);
				rr=(int)(sepsqrt*bin);		
				count[rr]=count[rr]+1.0;		
				rij=(RD[i]+RD[j])/2.0;
				srij=sepsqrt/rij;
				srr=(int)((srij-1)*bingof1);	
				
				if(srr < dim)
					g_of_1[srr]=g_of_1[srr]+1.0;
				if(i > j)
					rd_sum=rd_sum+pow(rij,(1.0*D));
			}
		}
	}
	count[0]=0.0;
}

void next2(int v[D])
{
        int ip;
        ip=0;
	while(v[ip] == 0)
	{
		ip++;
		if(ip == D)
			return;
	}

	while( ip <= D)
        {
                v[ip]=-v[ip];
                if(v[ip] < 0)
			break;	
		else if(v[ip] >= 0)
			ip=ip+1;
        }
}

void next(int v[D])
{
        int ip;
        ip=0;
        while(v[ip] <= Knum)
        {
                v[ip]=v[ip]+1;
                if(v[ip] <= Knum)
                        break;
                else if(v[ip] > Knum)
                {
                        v[ip]=0;
                        ip=ip+1;
                }
        }
}

int ipow(int base, int exp)
{
        int result = 1;
        while (exp)
        {
                if (exp & 1)
                        result *= base;
                exp >>= 1;
                base *= base;
        }
        return result;
}

void calculate_s_of_q()
{
	int i,j,q,u,k,v[D],count;
	double qq,qqsq,qqsqrt,qr;

	for(u=0;u<D;u++)
		v[u]=0;
	for(u=0;u<ipow((Knum+1),D);u++)
	{
		count=0;
		for(j=0;j<D;j++)
		{
			if(v[j] != 0)
				count++;
		}
		for(j=0;j<ipow(2,count);j++)
		{
			qqsq=0.0;
			for(k=0;k<D;k++)
				qqsq=qqsq+pow((1.0*v[k]),2.0);
			qqsqrt=sqrt(qqsq);
			q=(int)(qbin*qqsqrt);
			Sq_re[q]=0.0; Sqaa_re[q]=0.0;
			Sq_im[q]=0.0; Sqaa_im[q]=0.0;
			num[q]=num[q]+1.0; 
			numaa=0.0;
			for(i=0;i<N;i++)
			{	
				qr=0.0;
				for(k=0;k<D;k++)
					qr=qr+Kmin*v[k]*pos[k][i];
				Sq_re[q]=Sq_re[q]+cos(qr);
				Sq_im[q]=Sq_im[q]+sin(qr);
				if(RD[i] > Rbound)
				{
					Sqaa_re[q]=Sqaa_re[q]+cos(qr);
					Sqaa_im[q]=Sqaa_im[q]+sin(qr);
					numaa=numaa+1.0;
				}
			}
			Sq[q]=Sq[q]+pow(Sq_re[q],2.0)+pow(Sq_im[q],2.0);
			Sqaa[q]=Sqaa[q]+pow(Sqaa_re[q],2.0)+pow(Sqaa_im[q],2.0);
			Sqbb[q]=Sqbb[q]+pow((Sq_re[q]-Sqaa_re[q]),2.0)+pow((Sq_im[q]-Sqaa_im[q]),2.0);
			next2(v);
		}
		next(v);
	}
}

double vol_measure(double diam)
{
        double vol,g,dd,ag;
        dd=1.0*D;
        ag=D/2.0+1.0;
        vol=pow(pi,dd/2.0)*pow(diam,dd)/tgamma(ag);
        return vol;
}

double area_measure(double dist)
{
        double sarea,dd;
        int ag;
        dd=1.0*D;
        ag=D/2;
        sarea=2.0*pow(pi,dd/2.0)*pow(dist,(dd-1.0))/tgamma(ag);
        return sarea;
}

void vf_beyond_sphere()
{
	double locsq,loc,rp,l2;
	int i,k,rr;
	
	l2=1.0*RMAX_in*dr;
	for(i=0;i<SAMP;i++)
	{
		locsq=0.0;
		for(k=0;k<D;k++)
		{
			rp=boxl*(ran2(&seed)-0.50);
			locsq=locsq+pow(rp,2.0);
		}
		loc=sqrt(locsq);
		if(loc >= l2)
		{
			rr=(int)(loc*bin);	
			count_out[rr]=count_out[rr]+1.0;

		}
		if(loc < boxl2 && loc >= (boxl2-dr))
			count_in=count_in+1.0;
	}
}

void print_corr(int dev)
{
	FILE *fp;	
	int id;
	double tmp,nd,vol,vf_out,nd1;
	double sarea,pressure,sarea1;
	vol=pow(boxl,1.0*D);
	nd=1.0*N/vol;	
	sprintf(outfile1,"./correlation_out/real_sp_corr%dD_N%d_M%d_BS%1.4lfRmax%1.3lfRmin%1.3lfSWP%1.2lf_vf%1.4lfbf%1.5lfbin%d.dat",D,N,M,boxsize,rd,Rmin,swapprob,vf,buffer,bingof1);
	fp=fopen(outfile1,"w");
	fprintf(fp,"# Time = %d\n",(tstart+GAP2+(dev-1)*GAP2));
	sarea1=area_measure(boxl2-dr);
	for(id=0;id<RMAX;id++)
	{
		if(count[id] > 0)
		{
			tmp=1.0*id*dr;		
			sarea=area_measure(tmp);
			if(id < RMAX_in)
				fprintf(fp,"%lf\t%e\n",tmp,(double)(count[id]*bin)/sarea/nd/N/dev);
			else
			{
				if(count_out[id] > 0)
				{
					vf_out=sarea1*(count_out[id]/count_in);
					fprintf(fp,"%lf\t%e\n",tmp,(double)(count[id]*bin)/vf_out/nd/N/dev);
				}
			}
		}
	}
	fclose(fp);

	sprintf(outfile3,"./correlation_out/pressure%dD_N%d_M%d_BS%1.4lfRmax%1.3lfRmin%1.3lfSWP%1.2lf_vf%1.4lfbf%1.5lfbin%d.dat",D,N,M,boxsize,rd,Rmin,swapprob,vf,buffer,bingof1);
	fp=fopen(outfile3,"w");
	fprintf(fp,"# Time = %d NData=%d\n",(tstart+GAP2+(dev-1)*GAP2),dim-1);

	for(id=1;id<dim;id++)
	{
		tmp=1.0+1.0*id/bingof1;		
		sarea=area_measure(tmp);
	        nwg_of_1[id]=(double)(g_of_1[id]*bingof1)/sarea/dev/N/nd;
		fprintf(fp,"%lf\t%e\t%lf\n",tmp,nwg_of_1[id],rd_sum);
	}
	fclose(fp);

	sprintf(outfile2,"./correlation_out/S_k%dD_N%d_M%d_BS%1.4lfRmax%1.3lfRmin%1.3lfSWP%1.2lf_vf%1.4lf.dat",D,N,M,boxsize,rd,Rmin,swapprob,vf);
	fp=fopen(outfile2,"w");
	fprintf(fp,"# Time = %d\n",(tstart+GAP2+(dev-1)*GAP2));
	for(id=1;id<Kdim;id++)
	{
		if(Sq[id] > 0)
		{	
			tmp=1.0*id*Kmin/qbin;		
			fprintf(fp,"%lf\t%e\t%e\t%e\n",tmp,(double)(Sq[id])/num[id]/N,(double)(Sqaa[id])/num[id]/numaa,(double)(Sqbb[id])/num[id]/(N-numaa));
		}
	}
	fclose(fp);
}

void main()
{
	int k,t,tstrat;
	double tmpsq,tmp;
	tmpsq=0.0;
	scanf("%lf",&boxsize);
    scanf("%lf",&rd);
    scanf("%lf",&Rmin);
    scanf("%lf",&vf);
	scanf("%d",&nsteps);
	scanf("%lf",&Rskin);
	//boxsize=4.21743406; rd=1.20731; Rmin=0.8534; vf=0.2745; 
	//nsteps=4000000; Rskin=0.570000;
	GAP2=nsteps/1000;
	Rdif=rd-Rmin;
	Rbound=rd-Rdif*FRAC;
	boxl=boxsize; boxl2=((boxl)/2.0);
	cellsize=(rd+Rskin);
	Mcell=(int)(boxl/cellsize);
	tstart=nsteps-tnum*GAP2;
	dim=(buffer-1.0)*bingof1;
	for(k=0;k<D;k++)
		tmpsq=tmpsq+pow(boxl2,2.0);
	tmp=sqrt(tmpsq);
	Kmin=(2*pi/boxl);
	Kdim=(sqrt(1.0*D)*Knum*qbin)+1;

	RMAX=(int)(tmp*bin);
	RMAX_in=(int)(boxl2*bin);

	count=(double *) malloc ((RMAX+1) * sizeof(double));
	count_out=(double *) malloc ((RMAX+1) * sizeof(double));
	num=(double *) malloc (Kdim * sizeof(double));
	Sq_im=(double *) malloc (Kdim * sizeof(double));
	Sqaa_im=(double *) malloc (Kdim * sizeof(double));
	Sq_re=(double *) malloc (Kdim * sizeof(double));
	Sqaa_re=(double *) malloc (Kdim * sizeof(double));
	Sq=(double *) malloc (Kdim * sizeof(double));
	Sqaa=(double *) malloc (Kdim * sizeof(double));
	Sqbb=(double *) malloc (Kdim * sizeof(double));
	g_of_1 = (double *) malloc (dim * sizeof(double));
	nwg_of_1 = (double *) malloc (dim * sizeof(double));

	for(t=0;t<RMAX;t++)
	{
                count[t]=0.0;
		count_out[t]=0.0;
	}
	vf_beyond_sphere();

	for(t=0;t<tnum;t++)
	{
		read_input(t);
		calculate_g_of_r();
		calculate_s_of_q();
		print_corr(t+1);
	}
}

