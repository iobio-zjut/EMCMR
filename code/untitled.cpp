double MCSimulation::energyhbondnhoc2(point3f *decstr,int numseq)
{
	int i,j,k,l,m;
	point3d tp[20],ap[20],kp[20];
	point3d pd,pd2,pd3;
	double lamda;
	double diss[20];
	double totenergy=0;
	double totenergy2=0;

	static	double nhochhmean[22]={
		5.2172,8.7074,230.1397,4.2198,6.5391
		};
	static	double nhochhsigma[22]={
		0.3676,0.4375,10.2234,0.4080,0.3763
		};
	static	double nhocmeanval[][4]={//h-o h-o=c n-h-o n-h-o=c  
		2.85,89.0,110.5,199.5,//i+3
		2.00,147.0,159.0,160.0, //i+4  oldoh 
	//	2.83,89.0,110.0,201.5,//i+3 newoh
	//	2.00,148.0,159.0,155.0, //i+4

		2.00,155.0,164.0,180.0,//0
		2.00,155.0,164.0,180.0,//1
		2.00,151.0,163.0,192.0,//2
		2.00,151.0,163.0,192.0,//3

		};
	static	double nhocstdval[][4]={//h-o h-o=c n-h-o n-h-o=c 
		0.315504,7.697185,8.980366,7.932107,
		0.530859,10.582243,11.249764,25.360054,

		0.299730,11.770196,11.292558,68.955920,
		0.299730,11.770196,11.292558,68.955920,
		0.255088,12.376087,11.020081,69.165282,
		0.255088,12.376087,11.020081,69.165282,
		};

	double threshval[5]={19.113828,19.113828,20.723266,20.723266,16.118096};
	for(i=0;i<numseq;i++)
	{
		decstr[i].vpos=0;
		decstr[i].tpos=0;
		decstr[i].ssm='C';
		decstr[i].indl=-1;
		decstr[i].indr=-1;
		decstr[i].tpr=-1;
		decstr[i].tpl=-1;
	}
	for(i=1;i<numseq-4;i++)
	{
		j=i+3;	
		if((decstr[i].ss2=='E') || (decstr[j].ss2=='E' ))
			continue;
		l=0;
		for(k=0;k<4;k++)
		{
			if(decstr[i+k].ss2=='H') 
				l++;
		}
	/*	m=0;
		for(k=0;k<4;k++)
		{
			if(decstr[i+k].stype=='H') 
				m++;
		} */
		if(l==3 && m<2 && (decstr[i+1].ss2!='H' || decstr[i+2].ss2!='H'))
			continue;
		else if(l==2  && 
			((decstr[i+1].ss2=='H' && decstr[i+2].ss2=='H') || (decstr[i+1].ss2!='H' && decstr[i+2].ss2!='H') ))
			continue;
		if(l==0 ) continue;
		else if(l==1 ) continue;
		else if(l==2 ) lamda=0.05;
		else if(l==3) lamda=0.4;
		else if(l==4) lamda=1.5;
	//	else if(m==3) lamda=0.6;
	//	else if(m==4) lamda=1.0;
		else lamda=0.3;
		
		tp[1]=setv(decstr[i].x,decstr[i].y,decstr[i].z);
		ap[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
		tp[7]=minu(tp[1],ap[1]);
		diss[1]=norm(tp[7]); 
		if(diss[1]>6.6) continue;

		tp[0]=setv(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z);
		tp[2]=setv(decstr[i+1].x,decstr[i+1].y,decstr[i+1].z);
		ap[0]=setv(decstr[j-1].x,decstr[j-1].y,decstr[j-1].z);
		ap[2]=setv(decstr[j+1].x,decstr[j+1].y,decstr[j+1].z);
		kp[0]=minu(tp[2],ap[0]);
		diss[4]=norm(kp[0]);
		if(diss[4]<3.4 || diss[4]>4.2) continue;
		diss[0]=norm(minu(tp[0],ap[0]));
//		diss[1]=bf.norm(bf.minu(tp[0],ap[2]));			
//		diss[2]=bf.phi(tp[1].x,tp[1].y,tp[1].z,tp[2].x,tp[2].y,tp[2].z,ap[0].x,ap[0].y,ap[0].z,ap[1].x,ap[1].y,ap[1].z);
//		if(diss[2]>180) diss[2]-=180;
//		else diss[2]+=180;


//		tp[3]=bf.setv(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
//		tp[5]=bf.setv(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
//		ap[4]=bf.setv(decstr[j+1].ptn.x,decstr[j+1].ptn.y,decstr[j+1].ptn.z);
//		ap[5]=bf.setv(decstr[j].ptc.x,decstr[j].ptc.y,decstr[j].ptc.z);	
//		diss[3]=bf.norm(bf.minu(tp[3],ap[4]));
//		diss[4]=bf.norm(bf.minu(tp[5],ap[5]));

		//nhoc
		pd.x=decstr[j].pth.x-decstr[i].pto.x;
		pd.y=decstr[j].pth.y-decstr[i].pto.y;
		pd.z=decstr[j].pth.z-decstr[i].pto.z;
		diss[5]=norm(pd);
		if(diss[5]>=5.0 || diss[5]<1.6) continue;
 

		pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
		pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
		pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
		diss[6]=angv(pd,pd2)*degrad;
		if(diss[6]<70.0 || diss[6]>140.0) continue;

		pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
		pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
		pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
		diss[7]=angv(pd,pd3)*degrad;
				
		diss[8]=phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
			decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
		if(diss[8]<120.0 || diss[8]>280.0) continue;

		//second
		j=i+4;
		pd.x=decstr[j].pth.x-decstr[i].pto.x;
		pd.y=decstr[j].pth.y-decstr[i].pto.y;
		pd.z=decstr[j].pth.z-decstr[i].pto.z;
		diss[9]=norm(pd);
		if(diss[9]>=5.0 || diss[9]<1.6) continue;
 

		pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
		pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
		pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
		diss[10]=angv(pd,pd2)*degrad;
		if(diss[10]<100.0 || diss[10]>170.0) continue;

		pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
		pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
		pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
		diss[11]=angv(pd,pd3)*degrad;
		if(diss[11]<100.0) continue;
				
		diss[12]= phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
			decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
		if(diss[12]<100.0 || diss[12]>230.0) continue;

		double tval=0;
		for(k=0;k<1;k++)
		{		
			tval+=squgaussian(diss[k],nhochhsigma[k],nhochhmean[k]);
		}
		for(k=0;k<4;k++)
		{
			tval+=squgaussian(diss[5+k],nhocstdval[0][k],nhocmeanval[0][k]);
		}
		for(k=0;k<4;k++)
		{
			tval+=squgaussian(diss[9+k],nhocstdval[1][k],nhocmeanval[1][k]);
		}
		tval=exp(tval);
	//	printf("%3d %.3f\n",i,tval);
		if(tval>1e-7)
		{
			if(tval>1e-3) tval=1e-3;
			tval=lamda*(threshval[4]+log(tval));
			totenergy+=tval;
			for(k=0;k<=3;k++)
			{
				decstr[i+k].ssm='H';
			}
		}
	}
	
 
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////beta
	int inda,indb,indc,indd,inde,indf;
	int delseg;
	
	for(j=1;j<numseq-4;j++) 
	{
		for(k=j+3;k<numseq-1;k++)
		{	
			if( (decstr[j].ss2=='H' ) || (decstr[k].ss2=='H'))//diff
				continue;
			lamda=0.2;
			if((decstr[j].ss2=='E' ) && (decstr[k].ss2=='E' ))
				lamda+=0.1;
			else if((decstr[j].ss2=='E' ) || (decstr[k].ss2=='E' ))
				lamda+=0.2;
			else if(decstr[j].ss2=='E' && decstr[k].ss2=='E') 
				lamda+=1.6;
			else if(decstr[k].ss2=='E' || decstr[j].ss2=='E')
				lamda+=0.6;
			tp[1]=setv(decstr[j].x,decstr[j].y,decstr[j].z);
			ap[1]=setv(decstr[k].x,decstr[k].y,decstr[k].z);
			tp[7]=minu(tp[1],ap[1]);
			if(norm(tp[7])>8.0) continue;
		
			inda=j;indb=k;
			tp[4]=setv(decstr[inda].ptn.x,decstr[inda].ptn.y,decstr[inda].ptn.z);
			tp[5]=setv(decstr[inda+1].ptn.x,decstr[inda+1].ptn.y,decstr[inda+1].ptn.z);
			tp[9]=setv(decstr[inda-1].pto.x,decstr[inda-1].pto.y,decstr[inda-1].pto.z);
			tp[10]=setv(decstr[inda].pto.x,decstr[inda].pto.y,decstr[inda].pto.z);
			ap[4]=setv(decstr[indb].ptn.x,decstr[indb].ptn.y,decstr[indb].ptn.z);
			ap[5]=setv(decstr[indb+1].ptn.x,decstr[indb+1].ptn.y,decstr[indb+1].ptn.z);
			ap[9]=setv(decstr[indb-1].pto.x,decstr[indb-1].pto.y,decstr[indb-1].pto.z);
			ap[10]=setv(decstr[indb].pto.x,decstr[indb].pto.y,decstr[indb].pto.z);
			kp[0]=minu(tp[9],ap[4]);//oi-1 nj
			kp[1]=minu(tp[5],ap[10]);//ni+1 oj
			kp[2]=minu(tp[4],ap[9]);//ni oj-1
			kp[3]=minu(tp[10],ap[5]);//oi nj+1
			kp[4]=minu(tp[9],ap[5]);//oi-1 nj+1
			kp[5]=minu(tp[5],ap[9]);//ni+1 oj-1
			kp[6]=minu(tp[4],ap[10]);//ni oj
			kp[7]=minu(tp[10],ap[4]);//oi nj
			for(m=0;m<8;m++)
			{
				diss[m]=norm(kp[m]);
			}
			delseg=0;
			diss[9]=diss[0]+diss[1];
			diss[8]=diss[2]+diss[3];
			if(diss[8]<diss[9])
			{
				delseg=1;
				diss[9]=diss[8];
			}
			diss[8]=diss[4]+diss[5];
			if(diss[8]<diss[9])
			{
				delseg=2;
				diss[9]=diss[8];
			}
			diss[8]=diss[6]+diss[7];
			if(diss[8]<diss[9])
			{
				delseg=3;
				diss[9]=diss[8];
			}
			if(diss[2*delseg]>6.8) continue;
			if(diss[2*delseg+1]>6.8) continue;
			if(delseg==0)
			{
				indc=k;indd=j-1;inde=j+1;indf=k;
			}
			else if(delseg==1)
			{
				indc=j;indd=k-1;inde=k+1;indf=j;
			}
			else if(delseg==3)
			{
				indc=k;indd=j;inde=j;indf=k;
			}
			else if(delseg==2)
			{
				indc=k+1;indd=j-1;inde=j+1;indf=k-1;
			}

			pd.x=decstr[indc].pth.x-decstr[indd].pto.x;
			pd.y=decstr[indc].pth.y-decstr[indd].pto.y;
			pd.z=decstr[indc].pth.z-decstr[indd].pto.z;
			diss[0]=norm(pd);
			pd2.x=decstr[indd].ptc.x-decstr[indd].pto.x;
			pd2.y=decstr[indd].ptc.y-decstr[indd].pto.y;
			pd2.z=decstr[indd].ptc.z-decstr[indd].pto.z;
			diss[1]=angv(pd,pd2)*degrad;
			pd3.x=decstr[indc].pth.x-decstr[indc].ptn.x;
			pd3.y=decstr[indc].pth.y-decstr[indc].ptn.y;
			pd3.z=decstr[indc].pth.z-decstr[indc].ptn.z;
			diss[2]=angv(pd,pd3)*degrad;
	//		diss[3]=bf.phi(decstr[indc].ptn.x,decstr[indc].ptn.y,decstr[indc].ptn.z,decstr[indc].pth.x,decstr[indc].pth.y,decstr[indc].pth.z,
	//			decstr[indd].pto.x,decstr[indd].pto.y,decstr[indd].pto.z,decstr[indd].ptc.x,decstr[indd].ptc.y,decstr[indd].ptc.z);
			////////
			pd.x=decstr[inde].pth.x-decstr[indf].pto.x;
			pd.y=decstr[inde].pth.y-decstr[indf].pto.y;
			pd.z=decstr[inde].pth.z-decstr[indf].pto.z;
			diss[4]=norm(pd);				
			pd2.x=decstr[indf].ptc.x-decstr[indf].pto.x;
			pd2.y=decstr[indf].ptc.y-decstr[indf].pto.y;
			pd2.z=decstr[indf].ptc.z-decstr[indf].pto.z;
			diss[5]=angv(pd,pd2)*degrad;
			pd3.x=decstr[inde].pth.x-decstr[inde].ptn.x;
			pd3.y=decstr[inde].pth.y-decstr[inde].ptn.y;
			pd3.z=decstr[inde].pth.z-decstr[inde].ptn.z;
			diss[6]=angv(pd,pd3)*degrad;
	//		diss[7]=bf.phi(decstr[inde].ptn.x,decstr[inde].ptn.y,decstr[inde].ptn.z,decstr[inde].pth.x,decstr[inde].pth.y,decstr[inde].pth.z,
	//			decstr[indf].pto.x,decstr[indf].pto.y,decstr[indf].pto.z,decstr[indf].ptc.x,decstr[indf].ptc.y,decstr[indf].ptc.z);//n-h-o=c
	
			double tval=0;
			for(m=0;m<3;m++)
			{
				tval+=squgaussian(diss[m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
			}
			for(m=0;m<3;m++)
			{
				tval+=squgaussian(diss[4+m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
			}
			tval=exp(tval);

			if((delseg==0 && tval>5e-9)|| (delseg==1 && tval>5e-9)||(delseg==2 && tval>1e-9)|| (delseg==3 && tval>1e-9))
			{
				if(tval>1e-3) tval=1e-3;
				tval=lamda*(threshval[delseg]+log(tval))/threshval[delseg]*threshval[4];
				if(delseg==0 || delseg==2)//left of j
				{
					if(tval>decstr[j].vpos)
					{
						decstr[j].vpos=tval;
						if(decstr[j].indl>0) 
						{		
							if(decstr[decstr[j].indl].indr==j)//set to zero
							{
								decstr[decstr[j].indl].tpos=0;
								decstr[decstr[j].indl].indr=-1;
								if(decstr[decstr[j].indl].indl==-1)
									decstr[decstr[j].indl].ssm='C';
							}
						}
						decstr[j].indl=k;
						decstr[j].tpl=delseg;
						decstr[k].ssm='E';
						decstr[j].ssm='E';
					}
				}
				else if(delseg==1 || delseg==3)//right of j
				{
					if(tval>decstr[j].tpos)
					{
						decstr[j].tpos=tval;
						if(decstr[j].indr>0) 
						{
							if(decstr[decstr[j].indr].indl==j)//set to zero
							{
								decstr[decstr[j].indr].vpos=0;
								decstr[decstr[j].indr].indl=-1;
								if(decstr[decstr[j].indr].indr==-1)
									decstr[decstr[j].indr].ssm='C';
							}
						}
						decstr[j].indr=k;
						decstr[j].tpr=delseg;
						decstr[k].ssm='E';
						decstr[j].ssm='E';
					}
				}
				if(delseg==1 || delseg==2)//left of k
				{
					if(tval>decstr[k].vpos)
					{
						decstr[k].vpos=tval;
						if(decstr[k].indl>0) 
						{
							if(decstr[decstr[k].indl].indr==k)//set to zero
							{
								decstr[decstr[k].indl].tpos=0;
								decstr[decstr[k].indl].indr=-1;
								if(decstr[decstr[k].indl].indl==-1)
									decstr[decstr[k].indl].ssm='C';
							}
						}
						decstr[k].indl=j;
						decstr[k].tpl=delseg;
						decstr[j].ssm='E';
						decstr[k].ssm='E';
					}
				}
				else if(delseg==0 || delseg==3)//right of k
				{
					if(tval>decstr[k].tpos)
					{
						decstr[k].tpos=tval;
						if(decstr[k].indr>0) 
						{
							if(decstr[decstr[k].indr].indl==k)//set to zero
							{
								decstr[decstr[k].indr].vpos=0;
								decstr[decstr[k].indr].indl=-1;
								if(decstr[decstr[k].indr].indr==-1)
									decstr[decstr[k].indr].ssm='C';
							}
						}
						decstr[k].indr=j;
						decstr[k].tpr=delseg;
						decstr[j].ssm='E';
						decstr[k].ssm='E';
					}
				}	
			}
		}
	}
	double wthb=2.0;
	double wttot=3.0;
	if(protype==2) 
	{
		wthb=3.0;
		wttot=2.0;
	}
	for(i=0;i<numseq;i++)
	{
		totenergy2+=decstr[i].vpos;
		totenergy2+=decstr[i].tpos;
		if(decstr[i].indl!=-1 && decstr[i].indr!=-1)
		{
			totenergy2+=wthb;//beta 3.0 alpha 2.0
			if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==0 || decstr[i].tpr==1))
				totenergy2+=1.0;//+decstr[i].vpos+decstr[i].tpos;
			else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==2 || decstr[i].tpr==3))
				totenergy2+=1.5;//+decstr[i].vpos+decstr[i].tpos;
			else if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==2 || decstr[i].tpr==3))
				totenergy2-=0.0;
			else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==0 || decstr[i].tpr==1))
				totenergy2-=0.0;
		}
		else if(decstr[i].indl!=-1 || decstr[i].indr!=-1)
		{
			totenergy2+=1.5;
		}
	}
	return -(1.0*totenergy+wttot*totenergy2);//beta 1,2  alpha 1,3
}

bool tor2pos22(float xi,float yi,float zi,float xj,float yj,float zj,float xk,
		   float yk,float zk,float tang,float tleng,float tinner,float *xl,float *yl,float *zl)
{
	point3d p12,p23,e1,e2,e3;
	point3d p1,p2,p3;
	double tpangle,tcos,tsin,rmat[9];
	double q[4];
	p12.x=xi-xj;p12.y=yi-yj;p12.z=zi-zj;
	p23.x=xk-xj;p23.y=yk-yj;p23.z=zk-zj;
///*
	if(norm(p12)<epsilon*0.00001 || norm(p23)<epsilon*0.00001 || angv(p12,p23)<epsilon*0.001 || (PI-angv(p12,p23))<epsilon*0.001)
	{
		int imax=maxnormal(p23);
		if(imax==0)
		{
			yj+=0.01f;
		}
		else if(imax==1)
		{
			zj+=0.01f;
		}
		else
		{
			xj+=0.01f;
		}
		p12.x=xi-xj;p12.y=yi-yj;p12.z=zi-zj;
		p23.x=xk-xj;p23.y=yk-yj;p23.z=zk-zj;
		printf("make adjustment tor2pos22\n");
	}
//*/
	e2=unit(p23);
	e3=unit(p12);
	e1=prod(e2,e3);
	e1=unit(e1);
	if(norm(e1)<epsilon || norm(e2)<epsilon)
	{
		printf("wrong in tor2pos22 [%f %f %f] [%f %f %f] [%f %f %f]\n",xi,yi,zi,xj,yj,zj,xk,yk,zk);
		*xl=xk;*yl=yk;*zl=zk;
		return false;
	}
	p1=scal(e2,tleng);
	tpangle=(PI-tinner)/2.0;
	tcos=cos(tpangle);
	tsin=sin(tpangle);
	q[0]=tcos;q[1]=tsin*e1.x;q[2]=tsin*e1.y;q[3]=tsin*e1.z;
	q2rot(q,rmat);
	p2=mmat(rmat,p1);

	tpangle=tang/2.0;
	tcos=cos(tpangle);
	tsin=sin(tpangle);
	q[0]=tcos;q[1]=tsin*e2.x;q[2]=tsin*e2.y;q[3]=tsin*e2.z;
	q2rot(q,rmat);
	p3=mmat(rmat,p2);

	*xl=p3.x+xk;*yl=p3.y+yk;*zl=p3.z+zk;
	return true;
}

int maxnormal(point3d p1)
{
	int i;
	double tem,temp[3];	
	temp[0]=fabs(p1.x);
	temp[1]=fabs(p1.y);
	temp[2]=fabs(p1.z);
	i=0;
	tem=temp[0];
	if(temp[1]>tem)
	{
		i=1;
		tem=temp[1];
	}
	if(temp[2]>tem)
	{
		i=2;
		tem=temp[2];
	}
	return i;
}

void  q2rot(double *q,double rmax[])
{
	int i,j;
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			rmax[i*3+j]=0;
		}
	}
	rmax[0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
	rmax[1]=2*(q[1]*q[2]-q[0]*q[3]);
	rmax[2]=2*(q[1]*q[3]+q[0]*q[2]);
	rmax[3]=2*(q[1]*q[2]+q[0]*q[3]);
	rmax[4]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
	rmax[5]=2*(q[2]*q[3]-q[0]*q[1]);
	rmax[6]=2*(q[1]*q[3]-q[0]*q[2]);
	rmax[7]=2*(q[2]*q[3]+q[0]*q[1]);
	rmax[8]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
}
point3d mmat(double mat[3][3],point3d p1)
{
	point3d temp;
	temp.x=p1.x*mat[0][0]+p1.y*mat[0][1]+p1.z*mat[0][2];
	temp.y=p1.x*mat[1][0]+p1.y*mat[1][1]+p1.z*mat[1][2];
	temp.z=p1.x*mat[2][0]+p1.y*mat[2][1]+p1.z*mat[2][2];
	return temp;
}
point3d mmat(double mat[9],point3d p1)
{
	point3d temp;
	temp.x=p1.x*mat[0]+p1.y*mat[1]+p1.z*mat[2];
	temp.y=p1.x*mat[3]+p1.y*mat[4]+p1.z*mat[5];
	temp.z=p1.x*mat[6]+p1.y*mat[7]+p1.z*mat[8];
	return temp;
}

void ElectronDensity::select_points( vector<poseCoord > &pose, ObjexxFCL::FArray3D< float > densdata) {
//	ObjexxFCL::FArray3D< float > const & densdata = core::scoring::electron_density::getDensityMap().get_data();
	ObjexxFCL::FArray3D< std::complex<double> > Fdens, Frot;
	fourier::fft3(densdata, Fdens);

	// make rotationally averaged pose map
//	utility::vector1< core::Real > pose_1dspec;
	vector<float> pose_1dspec;
//	nRsteps_ =0 ;
	get_spectrum( pose, pose_1dspec);

//	if ( points_defined_ ) return; // points were predefined, don't change
	
	points_to_search_.clear();

	ObjexxFCL::FArray3D< double > rot;
	rot.dimension( densdata.u1(), densdata.u2(), densdata.u3() );
	map_from_spectrum( pose_1dspec, rot);
	fourier::fft3(rot, Frot);

	// FFT convolution
	int npts = densdata.u1()*densdata.u2()*densdata.u3();
	for ( int i=0; i< npts; ++i ) {
		Frot[i] = Fdens[i] * std::conj( Frot[i] );  // reuse
	}
	fourier::ifft3(Frot, rot);

	// mask asu
//	if ( symminfo_.enabled() ) {
//		symminfo_.mask_asu( rot , core::scoring::electron_density::getDensityMap(), 0);
//	}

	// sort points
//	utility::vector1< std::pair< numeric::xyzVector< core::Real >, core::Real > > point_score_pairs;
	int gridStep_=2;
	vector<std::pair< vector< float >, float > > point_score_pairs;
	for ( int z=1; z<=(int)densdata.u3(); z+=gridStep_ ) {
		for ( int y=1; y<=(int)densdata.u2(); y+=gridStep_ ) {
			for ( int x=1; x<=(int)densdata.u1(); x+=gridStep_ ) {
			//	numeric::xyzVector< core::Real > x_idx(x,y,z);
				vector<float> x_idx(3,0.0);
				x_idx[0] = x;
				x_idx[1] = y;
				x_idx[2] = z;
				float dens_value = -rot(x,y,z);
				std::pair< vector<float>, float> point_score_pair = std::make_pair(x_idx, dens_value);
				point_score_pairs.push_back(point_score_pair);
			}
		}
	}
	float minDistNative = 1e4;
	float point_radius_ = 3.0;
	int topNtrans_ = 5000;
	std::sort(point_score_pairs.begin(), point_score_pairs.end(), PointScoreComparator());
	for ( int i=0; i<point_score_pairs.size(); i++ ) {
		bool hasneighbor = false;
		vector<float> x_idx = point_score_pairs[i].first;
		vector<float> x_cart;
		x_cart[0] = x_idx[0];
		x_cart[1] = x_idx[1];
		x_cart[2] = x_idx[2];
		idx2cart( x_idx, x_cart );
		for ( int j=0; j< points_to_search_.size(); j++ ) {
			vector<float> x_idx_stored = points_to_search_[j];
			vector<float> x_cart_stored(3,0.0);
			x_cart_stored[0] = x_idx_stored[0];
			x_cart_stored[1] = x_idx_stored[1];
			x_cart_stored[2] = x_idx_stored[2];
//			(x_idx_stored[0],x_idx_stored[1],x_idx_stored[2]);
			idx2cart( x_idx_stored, x_cart_stored );
			float distance = sqrt(square_len(x_cart - x_cart_stored));
			if ( distance < point_radius_ ) hasneighbor = true;
		}
		if ( !hasneighbor ) {
			points_to_search_.push_back(point_score_pairs[i].first);
		//	if ( native_ ) {
		//		numeric::xyzVector< core::Real > x_cart;
		//		core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );
		//		core::Real distNative = symminfo_.min_symm_dist2(x_cart, native_com_);
		//		minDistNative = std::min( minDistNative, distNative );
		//	}
		}  
		if ( points_to_search_.size() >= topNtrans_ ) break;
	}


/*	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		core::scoring::electron_density::ElectronDensity mapdmp = core::scoring::electron_density::getDensityMap();
		mapdmp.set_data(rot);
		mapdmp.writeMRC( "filter.mrc" );
		//dump the points
		std::ofstream outpoints;
		outpoints.open("selectedpoints.txt", std::ofstream::app);
		for ( Size i=1; i<=points_to_search_.size(); i++ ) {
			numeric::xyzVector< core::Real > x_cart;
			numeric::xyzVector< core::Real > x_idx = points_to_search_[i];
			core::scoring::electron_density::getDensityMap().idx2cart( x_idx, x_cart );
			outpoints << "ATOM " << utility::to_string(i) << " " << x_cart[0] << " " << x_cart[1] << " " << x_cart[2] << std::endl;
		}
		outpoints.close();
	}  */

//	if ( native_ ) TR << "Closest point to native: " << std::sqrt(minDistNative) << std::endl;
}

void ElectronDensity::get_spectrum( vector<poseCoord > &pose, vector<float> &pose_1dspec) {
	vector<float> com(3,0.0);
	float extent = get_radius( pose, com );

	// grid spacing == delR
	delR_=2;
	int ngrid = int std::ceil( extent / delR_ + 2);
	vector<float>().swap(pose_1dspec);
	pose_1dspec=vector<float>(ngrid, 0.0);
	vector<float> pose_1dspec_for_nrsteps = pose_1dspec;

	float massSum=0.0;
/*	bool convolute_single_residue_ = false;
	// for fine point selection... generate pose_1dspec based on the middle residue of the pose.
	if ( convolute_single_residue_ == true ) {
	//	int midres = (pose.size()+1)/2;
		int midres= int((pose.size()/3 -1)/2) ;
	//	while ( pose.residue(midres).is_virtual_residue() ) {
	//		++midres;
	//	} 
	//	core::conformation::Residue const & rsd( pose.residue(midres) );
	//	core::conformation::Atom const & residue_CA( rsd.atom(2) );
		for ( int j=3*midres; j<= 3*midres+2; ++j ) {
	//		core::conformation::Atom const & atom( rsd.atom(j) );
			float binnum = ( pose[j].x_-residue_CA.xyz() ).length() / delR_ + 1.0;
			core::Real binnum = ( atom.xyz()-residue_CA.xyz() ).length() / delR_ + 1.0;
			core::Real fpart = binnum - std::floor(binnum);
			core::Size binint = (core::Size) std::floor(binnum);
			pose_1dspec[binint] += (1-fpart);
			pose_1dspec[binint+1] += (fpart);
		}
	} else { */
		// for coarse point selection... generate pose_1dspec based on whole pose
		for ( int i=0; i<pose.size(); ++i ) {
	//		core::conformation::Residue const & rsd( pose.residue(i) );
	//		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
	//		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
	//			core::conformation::Atom const & atom( rsd.atom(j) );
	//			core::Real binnum = ( atom.xyz()-com ).length() / delR_ + 1.0;
			vector<float> tmpC(3,0.0);
			tmpC[0] = com[0] - pose[i].x_[0];
			tmpC[1] = com[1] - pose[i].x_[1];
			tmpC[2] = com[2] - pose[i].x_[2];
			float binnum = sqrt(square_len( tmpC )) / delR_ + 1.0;
			float fpart = binnum - std::floor(binnum);
			int binint = (int) std::floor(binnum);
			pose_1dspec[binint] += (1-fpart);
			pose_1dspec[binint+1] += (fpart);
	//		}
		}
//	}

	// this is to calculate massSum via full pose for nRsteps_
	for ( int i=0; i<pose.size(); ++i ) {
//		core::conformation::Residue const & rsd( pose.residue(i) );
//		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
//		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
//			core::conformation::Atom const & atom( rsd.atom(j) );
			vector<float> tmpC(3,0.0);
			tmpC[0] = com[0] - pose[i].x_[0];
			tmpC[1] = com[1] - pose[i].x_[1];
			tmpC[2] = com[2] - pose[i].x_[2];
			float binnum = sqrt(square_len( tmpC))/ delR_ + 1.0;
			float fpart = binnum - std::floor(binnum);
			int binint = (int) std::floor(binnum);
			pose_1dspec_for_nrsteps[binint] += (1-fpart);
			pose_1dspec_for_nrsteps[binint+1] += (fpart);
			massSum += 1;
//		}
	}

	// now setting nRsteps_
	float  fracDens=0.7; // choose radius covering this fraction of density mass
	float running_total=0.0;
	for ( int i=0; i<(int)ngrid; ++i ) {
		running_total += pose_1dspec_for_nrsteps[i] / massSum;
		if ( running_total > fracDens ) {
			nRsteps_ = i;
			break;
		}
	}

	running_total=0.0;
	for ( int i=0; i<(int)ngrid; ++i ) {
		running_total += pose_1dspec[i] / massSum;
		pose_1dspec[i] /= (i*i); // normalize
//		  cout<< "spectrum " << i << ": " << pose_1dspec[i] << " " << pose_1dspec[i]*i*i/massSum << " " << running_total << std::endl;
	}
}


void ElectronDensity::map_from_spectrum( vector<float> &pose_1dspec, ObjexxFCL::FArray3D< double > &rot) {
	// ugly
	delR_=2;
	for ( int z=1; z<=(int)rot.u3(); ++z ) {
		for ( int y=1; y<=(int)rot.u2(); ++y ) {
			for ( int x=1; x<=(int)rot.u1(); ++x ) {
				vector<float> idxX(3,0.0), cartX(3,0.0);
				idxX[0] = x<=rot.u1()/2 ? x-1 : x-1-rot.u1();
				idxX[1] = y<=rot.u2()/2 ? y-1 : y-1-rot.u2();
				idxX[2] = z<=rot.u3()/2 ? z-1 : z-1-rot.u3();

				vector<float> fracX(3,0.0);
				fracX[0] = ( (float) idxX[0] ) / grid[0];
				fracX[1] = ( (float) idxX[1] ) / grid[1];
				fracX[2] = ( (float) idxX[2] ) / grid[2];
			//	cartX = f2c*fracX;
				MatrixTimesTransVector(f2c,fracX,cartX);
			//	idxoffset2cart( idxX, cartX );
				float d = sqrt(square_len(cartX)) / delR_;
				float fpart = d - std::floor(d);
				int dint = (int) std::floor(d) + 1;
				if ( dint<pose_1dspec.size() ) { // last entry is always 0 so this check is valid
					rot(x,y,z) = (1-fpart)*pose_1dspec[dint] + (fpart)*pose_1dspec[dint+1];
				} else {
					rot(x,y,z) = 0.0;
				}
			}
		}
	}

//	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
//		core::scoring::electron_density::ElectronDensity(rot,1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "spectrum.mrc" );
//	}
}

float ElectronDensity::get_radius( vector<poseCoord > & pose, vector<float> &com ){
	vector<float> centerCA(3,0.0);
	com = vector<float> (3,0.0);
	int nAtms=0;
	for ( int i=0; i< pose.size(); ++i ) {
	//	core::conformation::Residue const & rsd( pose.residue(i) );
	//	if ( rsd.aa() == core::chemical::aa_vrt ) continue;
	//	for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
	//		core::conformation::Atom const & atom( rsd.atom(j) );
			com[0] += pose[i].x_[0];
			com[1] += pose[i].x_[1];
			com[2] += pose[i].x_[2];
		//	if ( i==(pose.size()+1)/2  ) {
			if ( i== 3*int((pose.size()/3 -1)/2) +1) {
				int ii= 3*int((pose.size()/3 -1)/2) +1;
				centerCA[0] = pose[ii].x_[0];
				centerCA[1] = pose[ii].x_[1];
				centerCA[2] = pose[ii].x_[2];
			}
			nAtms++;
	//	}
	}
	com[0] /= nAtms;
	com[1] /= nAtms;
	com[2] /= nAtms;

//	if ( center_on_middle_ca_ ) {
//		com = centerCA;
//	}

	float maxSize = 0;
	for ( int i=0; i< pose.size(); ++i ) {
//		core::conformation::Residue const & rsd( pose.residue(i) );
//		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
//		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
//			core::conformation::Atom const & atom( rsd.atom(j) );
		vector<float> tmpCA(3,0.0);
		tmpCA[0] = com[0] - pose[i].x_[0];
		tmpCA[1] = com[1] - pose[i].x_[1];
		tmpCA[2] = com[2] - pose[i].x_[2];
		maxSize = std::max( maxSize, square_len(tmpCA));	
//			maxSize = std::max( maxSize, (com-atom.xyz()).length_squared());
//		}
	}
	maxSize = std::sqrt( maxSize ); // this the radius
	return maxSize;
}

void ElectronDensity::poseSphericalSamples(
	vector<poseCoord> &pose,
	ObjexxFCL::FArray3D< double > & sigR)
{
	using namespace core;

//	ElectronDensity &density = getDensityMap();

	int B=16;
	float delRsteps=2.0;
//	int nRsteps=nRsteps_;

	vector<float> reference_atm(3,0.0);
	vector<vector<float> > atmList;
	vector<float> all_K, all_C;

	vector<float> massSum(3,0.0), centerCA(3,0.0);

	// atom mask ... 3sigma from carbon
	float ATOM_MASK = 3.0 * sqrt( effectiveB / (2*M_PI*M_PI) );
//	float ATOM_MASK = 8.0;

	for ( int i=0; i< pose.size(); ++i ) {
//		conformation::Residue const & rsd( pose.residue(i) );
//		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
//		for ( Size j=1; j<= rsd.nheavyatoms(); ++j ) {
//			conformation::Atom const & atom( rsd.atom(j) );
//			atmList.push_back(atom.xyz());
		atmList.push_back(pose[i].x_);
		massSum[0] = massSum[0] + pose[i].x_[0];
		massSum[1] = massSum[1] + pose[i].x_[1];
		massSum[2] = massSum[2] + pose[i].x_[2];
//		massSum += atom.xyz();
		if ( i== 3*int((pose.size()/3 -1)/2) +1) {
			int ii= 3*int((pose.size()/3 -1)/2) +1;
			centerCA[0] = pose[ii].x_[0];
			centerCA[1] = pose[ii].x_[1];
			centerCA[2] = pose[ii].x_[2];
		}		

//			if ( i==(pose.size()+1)/2 && j==2 ) {
//				centerCA = atom.xyz();
//			}

//			chemical::AtomTypeSet const & atom_type_set( rsd.atom_type_set() );
//			std::string elt_i = atom_type_set[ rsd.atom_type_index( j ) ].element();
		std::string elt_i = pose[i].elt_;
		elt_i = elt_i[0];
		OneGaussianScattering sig_j = get_A( elt_i );

		float K_i = sig_j.k( effectiveB );
		all_K.push_back( K_i );
		all_C.push_back( sig_j.C( K_i ) );
//		}
	}
	int nAtms=atmList.size();
	massSum[0] /= nAtms; // center_of_mass = mass_sum / nAtms
	massSum[1] /= nAtms;
	massSum[2] /= nAtms;

//	if ( center_on_middle_ca_ ) {
//		massSum = centerCA;
//	}

	// precompute sines & cosines
	vector<float> cT,cG, sT,sG;
	cT.resize(2*B); cG.resize(2*B); sT.resize(2*B); sG.resize(2*B);
	for ( int t=0; t<2*B; ++t ) {
		float theta = (2.0*t-1.0)*M_PI/(4*B);
		sT[t] = sin(theta);
		cT[t] = cos(theta);
	}
	for ( int p=0; p<2*B; ++p ) {
		float phi = (2.0*p-2.0)*M_PI/(2*B);
		sG[p] = sin(phi);
		cG[p] = cos(phi);
	}

	//////////////////
	// pose -> spherical-sampled density
	// 1. one models each atom with a Gaussian sphere of density
	// 2. interpolate this calculated density in cencentric spherical shells
	// (extending out to D Ang in 1 Ang steps)
	//////////////////
//	if ( laplacian_offset_ != 0 ) {
//		TR << "Applying laplacian filter with offset of: " << laplacian_offset_ << " A" << std::endl;
//	}
	sigR.dimension( 2*B, 2*B, nRsteps_ );
	sigR = 0.0;

	// for each atom
	for ( int i=0; i<nAtms; ++i ) {
		float k=all_K[i];
		float C=all_C[i];

		atmList[i] -= massSum;

		float atomR = sqrt(square_len(atmList[i]));
		if ( atomR < 1e-5 ) {
			// uniform contribution to inner shells
			for ( int ridx=1; ridx<=nRsteps_; ++ridx ) {
				float atomD = ridx * delRsteps;
				if ( atomD < ATOM_MASK ) {
					float atomH = C * exp(-k*atomD*atomD); // <-- this is the place to calculate density
					for ( int t=1; t<=2*B; ++t ) {
						for ( int p=1; p<=2*B; ++p ) {
							sigR(p,t,ridx) += atomH;
						}
					}
				}
			}
			continue;
		}

		float beta = acos( atmList[i][2] / atomR );
		float gamma = atan2( atmList[i][0] , atmList[i][1] );   // x and y switched from usual convention

		float st1 = sin(beta);
		float sg1 = sin(gamma);
		float ct1 = cos(beta);
		float cg1 = cos(gamma);

		float laplacian_offset_ =0;
		if ( laplacian_offset_ != 0 ) {
			for ( int ridx=1; ridx<=nRsteps_; ++ridx ) {
				float shellR = ridx * delRsteps;
				for ( int t=1; t<=2*B; ++t ) {
					float minAtomD =  atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]+ct1*cT[t]);
					if ( minAtomD>ATOM_MASK*ATOM_MASK ) continue; // this just loops back to 2xB so we still get to sigR
					for ( int p=1; p<=2*B; ++p ) {
						float atomD = atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]*(sg1*sG[p]+cg1*cG[p])+ct1*cT[t]);
						if ( atomD < ATOM_MASK*ATOM_MASK ) {
							float atomH = C * exp(-k*atomD);
							sigR(p,t,ridx) += (-6 * atomH);

						}
					}
				}
			}
			// compute laplacian for surrounding coordinates
			for ( int xyz = 0; xyz < 3; ++xyz ) {
				for ( int lapl = 0; lapl < 2; ++lapl ) {
					reference_atm = atmList[i];
					atmList[i][xyz] = atmList[i][xyz] + ( ( (lapl==0) ? 1.0 : -1.0 ) * laplacian_offset_ );
					atomR = sqrt(square_len(atmList[i]));
					float beta = acos( atmList[i][2] / atomR );
					float gamma = atan2( atmList[i][0] , atmList[i][1] );   // x and y switched from usual convention
					float st1 = sin(beta);
					float sg1 = sin(gamma);
					float ct1 = cos(beta);
					float cg1 = cos(gamma);
					// residue index
					for ( int ridx=1; ridx<=nRsteps_; ++ridx ) {
						float shellR = ridx * delRsteps;
						for ( int t=1; t<=2*B; ++t ) {
							float minAtomD =  atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]+ct1*cT[t]);
							if ( minAtomD>ATOM_MASK*ATOM_MASK ) continue;
							for ( int p=1; p<=2*B; ++p ) {
								float atomD = atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]*(sg1*sG[p]+cg1*cG[p])+ct1*cT[t]);
								if ( atomD < ATOM_MASK*ATOM_MASK ) {
									float atomH = C * exp(-k*atomD);
									sigR(p,t,ridx) += atomH;
								}
							}
						}
					}
					// set atm back to original value
					atmList[i] = reference_atm;
				}
			}

		} else {
			for ( int ridx=1; ridx<=nRsteps_; ++ridx ) {
				float shellR = ridx * delRsteps;
				for ( Size t=1; t<=2*B; ++t ) {
					float minAtomD =  atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]+ct1*cT[t]);
					if ( minAtomD>ATOM_MASK*ATOM_MASK ) continue;
					for ( int p=1; p<=2*B; ++p ) {
						float atomD = atomR*atomR + shellR*shellR - 2*atomR*shellR*(st1*sT[t]*(sg1*sG[p]+cg1*cG[p])+ct1*cT[t]);
						if ( atomD < ATOM_MASK*ATOM_MASK ) {
							float atomH = C * exp(-k*atomD);
							sigR(p,t,ridx) += atomH;
						}
					}
				}
			}
		}
	} // loop through each atom

//	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
//		core::scoring::electron_density::ElectronDensity(sigR, 1.0, numeric::xyzVector< core::Real >(0,0,0), false).writeMRC( "Pose_sigR.mrc" );
//	}
}

void ElectronDensity::mapSphericalSamples (
	ObjexxFCL::FArray3D< double > &mapShellR,
	int nRsteps, float delR, int B,
	vector<float> center,
	float laplacian_offset
) {

	// make sure map is loaded
//	if ( !isLoaded ) {
//		TR.Fatal << "ElectronDensity::mapSHT called but no map is loaded!" << std::endl;
//		utility_exit();
//	}

//	numeric::xyzVector< core::Real > cartOffset, idxX, laplacian_cartOffset;
	vector<float> cartOffset(3,0.0), idxX(3,0.0), laplacian_cartOffset(3,0.0);
	float theta, phi, r;
	if ( coeffs_density_.u1()*coeffs_density_.u2()*coeffs_density_.u3() == 0 ) {
		spline_coeffs( density , coeffs_density_ );
	}

	mapShellR.dimension(2*B,2*B,nRsteps);

	// idx_com1: index coord
	vector<float> idx_com1 ( 3,center[0] );
	idx_com1[1] = center[1] ;
	idx_com1[2] = center[2] ;
	vector<float> idx_com2(3, center[0] + origin[0] - 1);
	idx_com2[1] = center[1] + origin[1] - 1 ;
	idx_com2[2] = center[2] + origin[2] - 1 ;
//	vector<float> idx_com2 ( center[0] + origin[0] - 1 ,
//		center[1] + origin[1] - 1 ,
//		center[2] + origin[2] - 1 );

	// resample
	// if laplacian
	if ( laplacian_offset != 0 ) {
		for ( int r_idx=1; r_idx<=nRsteps; ++r_idx ) {
			r = (float) delR*(r_idx);

			for ( int th_idx=1; th_idx<=2*B; ++th_idx ) {
				theta = (2.0*th_idx - 1.0) * M_PI / (4.0*B);

				for ( int phi_idx=1; phi_idx<=2*B; ++phi_idx ) {
					phi = (2.0*phi_idx - 2.0) * M_PI / (2.0*B);
					// reverse X/Y -- needed for spharm xform
					cartOffset[1] = r*cos(phi)*sin(theta);
					cartOffset[0] = r*sin(phi)*sin(theta);
					cartOffset[2] = r*cos(theta);
					laplacian_cartOffset = cartOffset;

					// apply laplacian filtering
					float ilaplsplinesum = 0;
					for ( int i = 0; i < 3; i++ ) {
						for ( int j = 0; j < 2; j++ ) {
							float tmp_lap = ( (j==0) ? 1.0 : -1.0);
							laplacian_cartOffset[i] = cartOffset[i] + ( tmp_lap * laplacian_offset );
							vector<float> fracX(3,0.0);
							MatrixTimesTransVector(c2f,laplacian_cartOffset,fracX);
						//	fracX= c2f*laplacian_cartOffset;
							idxX[0] = fracX[0]*grid[0] + idx_com1[0];
							idxX[1] = fracX[1]*grid[1] + idx_com1[1];
							idxX[2] = fracX[2]*grid[2] + idx_com1[2];
						//	idxX = numeric::xyzVector< core::Real >( fracX[0]*grid[0],
						//		fracX[1]*grid[1],
						//		fracX[2]*grid[2] );
						//	idxX = idxX + idx_com1;  // DAN!! is this right??
							ilaplsplinesum = ilaplsplinesum + interp_spline(coeffs_density_, idxX);
						}
					}

					// result of current point
					vector<float> fracX(3,0.0);
					MatrixTimesTransVector(c2f,cartOffset,fracX);
					idxX[0] = fracX[0]*grid[0] + idx_com1[0];
					idxX[1] = fracX[1]*grid[1] + idx_com1[1];
					idxX[2] = fracX[2]*grid[2] + idx_com1[2];
				//	idxX = numeric::xyzVector< core::Real >( fracX[0]*grid[0],
				//		fracX[1]*grid[1],
				//		fracX[2]*grid[2] );
				//	idxX = idxX + idx_com1;

					// laplacian is -6 * current_point + sum of the +-1 of surrounding axes
					mapShellR(phi_idx, th_idx, r_idx) = (-6 * interp_spline( coeffs_density_ , idxX ) ) + ilaplsplinesum;
				}
			}
		}
	} else {
		for ( int r_idx=1; r_idx<=nRsteps; ++r_idx ) {
			r = (float) delR*(r_idx);

			for ( int th_idx=1; th_idx<=2*B; ++th_idx ) {
				theta = (2.0*th_idx - 1.0) * M_PI / (4.0*B);

				for ( int phi_idx=1; phi_idx<=2*B; ++phi_idx ) {
					phi = (2.0*phi_idx - 2.0) * M_PI / (2.0*B);
					// reverse X/Y -- needed for spharm xform
					cartOffset[1] = r*cos(phi)*sin(theta);
					cartOffset[0] = r*sin(phi)*sin(theta);
					cartOffset[2] = r*cos(theta);

					vector<float> fracX(3,0.0);
					MatrixTimesTransVector(c2f,cartOffset,fracX);
				//	 = c2f*cartOffset;
					idxX[0] = fracX[0]*grid[0] + idx_com1[0];
					idxX[1] = fracX[1]*grid[1] + idx_com1[1];
					idxX[2] = fracX[2]*grid[2] + idx_com1[2];					
				//	idxX = numeric::xyzVector<core::Real>( fracX[0]*grid[0],
				//		fracX[1]*grid[1],
				//		fracX[2]*grid[2] );

				//	idxX = idxX + idx_com1;

					mapShellR(phi_idx, th_idx, r_idx) = interp_spline( coeffs_density_ , idxX );
				}
			}
		}
	}
}

// do the main search over the map
void ElectronDensity::density_grid_search (
	int pose_idx,
	vector<poseCoord> & pose,
	RBfitResultDB & results
) {
//	core::scoring::electron_density::ElectronDensity &density = core::scoring::electron_density::getDensityMap();
//	numeric::xyzVector<int> grid = density.getGrid();

	nRsteps_=0;
	delR_=2;
	B_ = 16;	
	select_points(pose, density);
	// allocate space for SHT
	fourier::SHT SOFT(B_, nRsteps_);

	// get com of pose
	vector<vector<float>> rot;
	vector<float> pretrans, posttrans;
	get_radius( pose, pretrans );

	pretrans[0]=-1.0*pretrans[0];
	pretrans[1]=-1.0*pretrans[1];
	pretrans[2]=-1.0*pretrans[2];

//	runtime_assert( points_to_search_.size() >= 1 ); // sanity check

	// get pose SPHARM
	ObjexxFCL::FArray3D< double > poseSig, poseCoefR, poseCoefI;
	poseSphericalSamples( pose, poseSig);
	SOFT.sharm_transform( poseSig, poseCoefR, poseCoefI );
	SOFT.sph_standardize( poseCoefR, poseCoefI );

	ObjexxFCL::FArray3D< double > mapSig, mapCoefR, mapCoefI;
	for ( int i=0; i< points_to_search_.size(); ++i ) {
		// get cartesian coords of ths point
		density.idx2cart( points_to_search_[i], posttrans );

		density.mapSphericalSamples( mapSig, nRsteps_, delR_, B_, points_to_search_[i], laplacian_offset_ );
		SOFT.sharm_transform( mapSig, mapCoefR, mapCoefI );
		SOFT.sph_standardize( mapCoefR, mapCoefI );

		// get correlation
		ObjexxFCL::FArray3D< double > so3_correlation;
		SOFT.so3_correlate(so3_correlation, mapCoefR,mapCoefI,  poseCoefR,poseCoefI);


		//core::Size nperRot = 25*std::max((core::Size)3, topNfilter_ / points_to_search_.size());
		//core::Size nperRotCl = std::max((core::Size)3, topNfilter_ / points_to_search_.size());

		// we initially oversample since the set is so clustered
		// no matter how many we want, don't take more than 1/8 of everything (which still might be a lot)
		core::Size nperRot = std::min( 100*max_rot_per_trans_ , B_*B_*B_);
		core::Size nperRotCl = max_rot_per_trans_;

		RBfitResultDB local_results( nperRot );

		if ( normscores_ ) {
			core::Real correl_sum=0.0, correl_sum2=0.0;
			for ( core::Size j=0; j<8*B_*B_*B_; ++j ) {
				correl_sum += so3_correlation[j];
				correl_sum2 += so3_correlation[j]*so3_correlation[j];
			}
			correl_sum /= 8*B_*B_*B_;
			correl_sum2 = std::sqrt( correl_sum2/(8*B_*B_*B_) - correl_sum*correl_sum );
			for ( core::Size j=0; j<8*B_*B_*B_; ++j ) {
				so3_correlation[j] = (so3_correlation[j] - correl_sum) / correl_sum2;
			}
		}

		for ( core::Size j=0; j<8*B_*B_*B_; ++j ) {
			if ( local_results.to_add_element( so3_correlation[j] ) ) {
				SOFT.idx_to_rot(j , rot);
				local_results.add_element( RBfitResult( pose_idx, so3_correlation[j], rot, pretrans, posttrans ) );
			}
		}

		do_filter( local_results );
		while ( local_results.size() > nperRotCl ) local_results.pop();

		core::Size nclust = local_results.size();

		core::Real bestrms=1e6,bestscore=0;
		while ( local_results.size() > 0 ) {
			RBfitResult sol_i = local_results.pop();

			if ( native_ ) {
				core::pose::PoseOP posecopy ( new core::pose::Pose(pose) );
				apply_transform( *posecopy, sol_i );
				core::pose::addVirtualResAsRoot( *posecopy );
				core::Real rms_i = get_rms(native_, posecopy, symminfo_);
				if ( rms_i < bestrms ) {
					bestrms = rms_i; bestscore = sol_i.score_;
				}
			}

			results.add_element( sol_i );
		}

		core::Real minDistNative=1e6;
		if ( native_ ) {
			core::Real distNative = (posttrans-native_com_).length_squared();
			minDistNative = std::min( minDistNative, distNative );
		}

		if ( std::sqrt(minDistNative) <5.0 ) {
			TR << "[" << i << "/" << points_to_search_.size() << "]" << " nmdls " << nclust << " pointrms " << std::sqrt(minDistNative)
				<< " rms " << bestrms << " score " << bestscore << std::endl;
		}
		if ( i%100 == 0 ) {
			TR << "[" << i << "/" << points_to_search_.size() << "] " << results.top().score_ << " / " << results.size() << std::endl;
		}
	}

	//TR << "[" << points_to_search_.size() << "/" << points_to_search_.size() << "]" << std::endl;
}

bool mcfragsweepaaa(vector<point3f> &tmstr,int numseq,int numshelix)
{
	int i;
	int trandpos;
	double trand2,tlength;
	int leng1,leng2,leng3,indbin;
	double rmat[9];
	point3d tp[12];
	double tnewenergy;

    vector<sssegment> sse;
    int numsse=0;
    vector<int> alphasind;
    vector<int> alphaind;
    vector<int> betaind;
    genesse(tmstr,numseq,sse,numsse);
//    genessex(tmstr,numseq,sse,numsse,alphasind);
    calcabind(alphasind,alphaind,betaind,numsee,sse);	
//	point3f *tmstr=new point3f[numseq];
//	BasicFunc bf;
//	ParseSeq ps;
	bool flagclash=false;
	bool flagtor;
	double threshclash=0.01;
//	int numclash;
/*    int tothelix=0;
    for(i=0;i<numsse;i++)
    {
        if(sse[i].ss=='H') tothelix++;
    }
//    cout<<"CC: "<<tothelix<<endl;
    if(tothelix==0)
    {
//        nummov[13][1]++;
        return false;
    }	*/
    int numsalpha = alphasind.size();
	if(numsalpha==0)
	{
//		summcaaa[2]++;
//		delete[]tmstr;
		return false;
	}
do
{
//	for(i=0;i<numseq;i++)
//	{
//		tmstr[i]=decstr[i];
// 	}
//	memcpy(tmstr,decstr,numseq*sizeof(point3f));
//	accdispos[0]=dispos[0];
//	for(i=1;i<numseq;i++)
//	{
//		accdispos[i]=accdispos[i-1]+dispos[i];
//	}	
		//only in coil
//		trandpos=findindpos(accdispos,accdispos[numseq-1]*rand()/double(RAND_MAX+1.0),numseq);
//		if(tmstr[trandpos].stype!='C' || trandpos==numseq-1)
//		{
//			continue;
//		}

//		trandpos=int((numsalpha)*(rand()/double(RAND_MAX+1.0)));
//		trandpos=int((numsalpha)*Random());
		trandpos=numshelix;
		leng1=sse[alphasind[trandpos]].term-sse[alphasind[trandpos]].init+1;
		leng2=sse[alphasind[trandpos]+1].term-sse[alphasind[trandpos]+1].init+1;
		leng3=sse[alphasind[trandpos]+2].term-sse[alphasind[trandpos]+2].init+1;
		if(leng1>4 && leng3>4 && leng3<20)
		{
			indbin=0;
			trand2=randf0and1();
			for(i=0;i<30;i++)
			{
				if(trand2<=angleaa[i])
				{
					indbin=i;
					break;
				}
			}
			trand2=(indbin+randf0and1())/30.0*PI;
		}
		else 
			trand2=randf0and1()*PI;
		tp[0]=setv(tmstr[sse[alphasind[trandpos]].init].x,tmstr[sse[alphasind[trandpos]].init].y,
			tmstr[sse[alphasind[trandpos]].init].z);
		tp[1]=setv(tmstr[sse[alphasind[trandpos]].term].x,tmstr[sse[alphasind[trandpos]].term].y,
			tmstr[sse[alphasind[trandpos]].term].z);
		tp[2]=setv(tmstr[sse[alphasind[trandpos]+2].init].x,tmstr[sse[alphasind[trandpos]+2].init].y,
			tmstr[sse[alphasind[trandpos]+2].init].z);//for new begin more flexible
		tp[3]=setv(tmstr[sse[alphasind[trandpos]+2].term].x,tmstr[sse[alphasind[trandpos]+2].term].y,
			tmstr[sse[alphasind[trandpos]+2].term].z);
		tp[4]=minu(tp[3],tp[2]);
		tlength=norm(tp[4]);
		tp[5]=rana(trand2);
		tp[6]=setv(0,0,1);
		tp[7]=minu(tp[0],tp[1]);
		tp[8]=unit(tp[7]);
		v2rot(tp[8],rmat);
		tp[9]=rotv(tp[5],rmat);
		tp[10]=scalx(tp[9],tlength);
		tp[11]=addv(tp[2],tp[10]);//for new end
		mcfragsweepCCD4(tmstr,numseq,sse[alphasind[trandpos]+1].init,sse[alphasind[trandpos]+1].term,
			tp[2],sse[alphasind[trandpos]+2].init,tp[11],sse[alphasind[trandpos]+2].term);
		flagtor=tor2str(tmstr,numseq,3);
		if(!flagtor)
		{
			printf("tor2str wrong in aaa\n");
		}
//	numclash=ps.calcaclash(tmstr,numseq);
//	if(numclash<numseq*threshclash)
//	{
//		flagclash=true;
//	}
//	else flagclash=false;
		flagclash=true;
}while(!flagclash);
	
//	int numclash;
//	numclash=ps.calcaclash(tmstr,numseq);
//	if(numclash>numseq*threshclash)
//	{
//		summcome[2]++;
//		delete[]tmstr;
//		return false;
//	}
	//////////make decision
//	bool flagev=evcheck(tmstr,numseq);
//	if(!flagev) numev3[10][3]++;
//	if(flagev)
//	{
//		summcaaa[2]++;
//		return false;
//	}
	return false;
}

bool mcfragsweepCCD4(vector<point3f> &tmstr2,int numseq,int lps,int lpe,point3d pt1,int ind1,point3d pt2,int ind2)
{
//	int i;
//	BasicFunc bf;
//	ParseSeq ps;
	int trandpos;
//	int indphi,indpsi;
	bool flagphi,flagtor;
	int flagpt,indp[2];
	indp[0]=ind1;indp[1]=ind2;
	point3d tp[2]; 
	tp[0].x=pt1.x;tp[0].y=pt1.y;tp[0].z=pt1.z;
	tp[1].x=pt2.x;tp[1].y=pt2.y;tp[1].z=pt2.z;
//	point3f *tmstr=new point3f[numseq];
	point3d pcur,p12,p13;
	int numiter=0;
	double tdist[2],ttheta,tphi,tpsi,tinner;
//	for(i=0;i<numseq;i++)
//	{
//		tmstr[i]=decstr[i];
//	}
//	memcpy(tmstr2,decstr,numseq*sizeof(point3f));
	flagpt=0;
//	trandpos=lps;
	int threshiter=25*(lpe-lps+1);
	if(threshiter>600) threshiter=600;
	do{
		tdist[0]=10000;
		tdist[1]=10000;
//		trandpos=lps+int((lpe-lps+1)*(rand()/double(RAND_MAX+1.0)));
		trandpos=lps+int((lpe-lps+1)*randf0and1());
//		if(trandpos>lpe)
//		{
//			printf("ccd3 pos is large %d %d %d\n",lps,lpe,trandpos);
//			trandpos=lpe;
//		}
//		trandpos=lps+(trandpos+1-lps)%(lpe-lps+1);
	//	if(genrand()%2==0)
		if(randf0and1()<0.5)
		{
			flagphi=true;//change phi n ca
		}
		else
		{
			flagphi=false;//change psi ca c in pos+1
		}
		if(tmstr2[trandpos].ss2=='H' && tmstr2[trandpos+1].ss2=='H')
		{
			numiter++;
			continue;
		}
		else if(tmstr2[trandpos].ss2=='H')
		{
			flagphi=false;
		}
		else if(tmstr2[trandpos].ss2=='E' && tmstr2[trandpos+1].ss2=='E')
		{
			numiter++;
			continue;
		}
		else if(tmstr2[trandpos].ss2=='E')
		{
			flagphi=false;
		}
		flagpt=(flagpt+1)%2;
		if(flagphi)
		{
			p12=setv(tmstr2[trandpos].x-tmstr2[trandpos].ptn.x,tmstr2[trandpos].y-tmstr2[trandpos].ptn.y,
				tmstr2[trandpos].z-tmstr2[trandpos].ptn.z);
			p13=setv(tmstr2[indp[flagpt]].x-tmstr2[trandpos].ptn.x,tmstr2[indp[flagpt]].y-tmstr2[trandpos].ptn.y,
				tmstr2[indp[flagpt]].z-tmstr2[trandpos].ptn.z);
			tinner=angv(p12,p13);
			if(tinner<epsilon || PI-tinner<epsilon || norm(p12)<epsilon)//in one line no affect when rotate
			{
				numiter++;
			//	printf("oneline 1\n");
				continue;
			}
			pcur.x=tmstr2[indp[flagpt]].x;pcur.y=tmstr2[indp[flagpt]].y;pcur.z=tmstr2[indp[flagpt]].z;
			p13.x=tmstr2[trandpos].x;p13.y=tmstr2[trandpos].y;p13.z=tmstr2[trandpos].z;
			p12.x=tmstr2[trandpos].ptn.x;p12.y=tmstr2[trandpos].ptn.y;p12.z=tmstr2[trandpos].ptn.z;
			ttheta=phi(pcur.x,pcur.y,pcur.z,p12.x,p12.y,p12.z,p13.x,p13.y,p13.z,tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
			tphi=tmstr2[trandpos].tor[2]+ttheta;
			if(tphi>360.0)
			{
				tphi-=360.0;
			}
//			tpsi=tmstr2[trandpos+1].tor[0];
//			indphi=int(tphi/18.0);
//			indpsi=int(tpsi/18.0);
//			if(!tphipsi[indpsi][indphi])
//			{
//			//	continue;
//			}
			tmstr2[trandpos].tor[2]=tphi;	
		}
		else
		{
			p12=setv(tmstr2[trandpos].ptc.x-tmstr2[trandpos].x,tmstr2[trandpos].ptc.y-tmstr2[trandpos].y,
				tmstr2[trandpos].ptc.z-tmstr2[trandpos].z);
			p13=setv(tmstr2[indp[flagpt]].x-tmstr2[trandpos].x,tmstr2[indp[flagpt]].y-tmstr2[trandpos].y,
				tmstr2[indp[flagpt]].z-tmstr2[trandpos].z);
			tinner=angv(p12,p13);
			if(tinner<epsilon || PI-tinner<epsilon || norm(p12)<epsilon)//in one line no affect when rotate
			{
				numiter++;
			//	printf("oneline 2\n");
				continue;
			}
			pcur.x=tmstr2[indp[flagpt]].x;pcur.y=tmstr2[indp[flagpt]].y;pcur.z=tmstr2[indp[flagpt]].z;
			p12.x=tmstr2[trandpos].x;p12.y=tmstr2[trandpos].y;p12.z=tmstr2[trandpos].z;
			p13.x=tmstr2[trandpos].ptc.x;p13.y=tmstr2[trandpos].ptc.y;p13.z=tmstr2[trandpos].ptc.z;
			ttheta=phi(pcur.x,pcur.y,pcur.z,p12.x,p12.y,p12.z,p13.x,p13.y,p13.z,tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
			tpsi=tmstr2[trandpos+1].tor[0]+ttheta;
			if(tpsi>360.0)
			{
				tpsi-=360.0;
			}
//			tphi=tmstr2[trandpos].tor[2];
//			indphi=int(tphi/18.0);
//			indpsi=int(tpsi/18.0);
//			if(!tphipsi[indpsi][indphi])
//			{
//			//	continue;//go to the endpos to check tdist and numiter
//			}
			tmstr2[trandpos+1].tor[0]=tpsi;
		}		
		numiter++;
		if(trandpos==0)
			flagtor=ps.tor2str(tmstr2,numseq,3);
		else
			flagtor=ps.tor2strp(tmstr2,numseq,trandpos);
		if(!flagtor)
		{
			printf("tor2str wrong in CCD4 %d\n",trandpos);
		}
	pcur.x=tmstr2[ind1].x-tp[0].x;
	pcur.y=tmstr2[ind1].y-tp[0].y;
	pcur.z=tmstr2[ind1].z-tp[0].z;
	tdist[0]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
	pcur.x=tmstr2[ind2].x-tp[1].x;
	pcur.y=tmstr2[ind2].y-tp[1].y;
	pcur.z=tmstr2[ind2].z-tp[1].z;
	tdist[1]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
//	printf("%d %d %f %f %f\n",numiter,threshiter,tdist[0],tdist[1],tdist[2]);
	}while(numiter<threshiter && (tdist[0]>6.50 || tdist[1]>2.50));
//	printf("%d [%d %d %d] pos %d %d dist %f theta %f\n",numiter,lps,lpe,lpt,trandpos, flagphi, tdist,ttheta);
//	*newenergy=calcrmsdenergy(tmstr2,numseq);
/*	if(numiter<threshiter)
	{

//
//		for(i=0;i<numseq;i++)
//		{
//			decstr[i]=tmstr[i];
//		}
		memcpy(decstr,tmstr2,numseq*sizeof(point3f));
//		delete[]tmstr;
		return true;
	}
	else
	{
//		for(i=0;i<numseq;i++)
//		{
//			decstr[i]=tmstr[i];
// 		}
		memcpy(decstr,tmstr2,numseq*sizeof(point3f));
//		delete[]tmstr;
		return false;
	}  */
	return true;
}

void calcabind(vector<int> &alphasind, vector<int> &alphaind, vector<int> &betaind,int numsee,vector<sssegment> &sse)
{
	int i;
/*	if(alphaind)
	{
		delete[]alphaind;
		alphaind=NULL;
	}
	if(betaind)
	{
		delete[]betaind;
		betaind=NULL;
	}
	if(alphasind)
	{
		delete[]alphasind;
		alphasind=NULL;
	} */
	int numsalpha=0;
	int numalpha=0;
	int numbeta=0;
	alphaind = vector<int>(numsse,0);
	betaind = vector<int>(numsse,0);
	alphasind = vector<int>(numsse,0);
//	alphaind=new int[numsse];
//	alphasind=new int[numsse];
//	betaind=new int[numsse];
	for(i=0;i<numsse;i++)
	{
		if(sse[i].ss=='H')
			alphaind[numalpha++]=i;
		else if(sse[i].ss=='E')
			betaind[numbeta++]=i;
	}
	for(i=0;i<numalpha-1;i++)
	{
		if(alphaind[i+1]==alphaind[i]+2)
		{
			alphasind[numsalpha++]=alphaind[i];
		}
	}
}

bool tor2str(vector<point3f> &decstr,int seqnum,int type)
{
	int i;
	point3s pt,pn,pc;
	BasicFunc bf;
	bool flagts;
	bool flagwhole=true;
if(type==1)//all 0 n 1 ca 2 c
{
	if(decstr[1].leng<0)
	{
		decstr[1].leng=float(lennca);
	}
	if(decstr[2].leng<0)
	{
		decstr[2].leng=float(lencac);
	}
	if(decstr[2].angl<0)
	{
		decstr[2].angl=float(angncac);
	}
	decstr[0].x=0;decstr[0].y=0;decstr[0].z=0;
	decstr[1].x=decstr[1].leng;decstr[1].y=0;decstr[1].z=0;
	decstr[2].x=decstr[1].leng-decstr[2].leng*cos(decstr[2].angl*raddeg);
	decstr[2].y=decstr[2].leng*sin(decstr[2].angl*raddeg);decstr[2].z=0;
	for(i=3;i<seqnum;i++)
	{
		flagts=tor2pos22(decstr[i-3].x,decstr[i-3].y,decstr[i-3].z,decstr[i-2].x,decstr[i-2].y,decstr[i-2].z,
			decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i].phi*raddeg,decstr[i].leng,decstr[i].angl*raddeg,&pt.x,&pt.y,&pt.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinates %d\n",i);
		}
		decstr[i].x=pt.x;decstr[i].y=pt.y;decstr[i].z=pt.z;
	}
}
else if(type==6)//ca 
{
	if(decstr[1].leng<0)
	{
		decstr[1].leng=float(lencaca);
	}
	if(decstr[2].leng<0)
	{
		decstr[2].leng=float(lencaca);
	}
	if(decstr[2].angl<0)
	{
		decstr[2].angl=106.422f;
	}
	decstr[0].x=0;decstr[0].y=0;decstr[0].z=0;
	decstr[1].x=decstr[1].leng;decstr[1].y=0;decstr[1].z=0;
	decstr[2].x=decstr[1].leng-decstr[2].leng*cos(decstr[2].angl*raddeg);
	decstr[2].y=decstr[2].leng*sin(decstr[2].angl*raddeg);decstr[2].z=0;
	for(i=3;i<seqnum;i++)
	{
		flagts=tor2pos22(decstr[i-3].x,decstr[i-3].y,decstr[i-3].z,decstr[i-2].x,decstr[i-2].y,decstr[i-2].z,
			decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i].phi*raddeg,decstr[i].leng,decstr[i].angl*raddeg,&pt.x,&pt.y,&pt.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinates %d\n",i);
		}
		decstr[i].x=pt.x;decstr[i].y=pt.y;decstr[i].z=pt.z;
	}
}
else if(type==3)
{
	if(decstr[0].len[1]<0)
	{
		decstr[0].len[1]=float(lennca);
	}
	if(decstr[0].len[2]<0)
	{
		decstr[0].len[2]=float(lencac);
	}
	if(decstr[0].ang[2]<0)
	{
		decstr[0].ang[2]=float(angncac);
	}
	decstr[0].ptn.x=0;decstr[0].ptn.y=0;decstr[0].ptn.z=0;
	decstr[0].x=decstr[0].len[1];decstr[0].y=0;decstr[0].z=0;
	decstr[0].ptc.x=decstr[0].len[1]-decstr[0].len[2]*cos(decstr[0].ang[2]*raddeg);
	decstr[0].ptc.y=decstr[0].len[2]*sin(decstr[0].ang[2]*raddeg);decstr[0].ptc.z=0;
	if(decstr[0].tor[0]>=0 && decstr[0].tor[2]>=0)
	{
		flagts=tor2pos22(0,0,0,lennca,0,0,lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
			decstr[0].tor[0]*raddeg,float(lencn),float(angcacn)*raddeg,&pn.x,&pn.y,&pn.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinates n %d\n",0);
		}
		decstr[0].ptn.x=pn.x;decstr[0].ptn.y=pn.y;decstr[0].ptn.z=pn.z;
		if(decstr[0].tor[1]<0) 
			decstr[0].tor[1]=180;
		flagts=bf.tor2pos22(lennca,0,0,lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
			decstr[0].ptn.x,decstr[0].ptn.y,decstr[0].ptn.z,
			decstr[0].tor[1]*raddeg,float(lennca),angcnca*raddeg,&pt.x,&pt.y,&pt.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinates ca %d\n",0);
		}
		decstr[0].x=pt.x;decstr[0].y=pt.y;decstr[0].z=pt.z;
		flagts=bf.tor2pos22(lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
			decstr[0].ptn.x,decstr[0].ptn.y,decstr[0].ptn.z,decstr[0].x,decstr[0].y,decstr[0].z,
			decstr[0].tor[2]*raddeg,float(lencac),angncac*raddeg,&pc.x,&pc.y,&pc.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinates c %d\n",0);
		}
		decstr[0].ptc.x=pc.x;decstr[0].ptc.y=pc.y;decstr[0].ptc.z=pc.z;
	}
	for(i=1;i<seqnum;i++)
	{
		if(decstr[i].tor[0]<0) 
		{
			if(decstr[i].stype=='E')
				decstr[i].tor[0]=120;
			else if(decstr[i].stype=='H')
				decstr[i].tor[0]=300;
			else
				decstr[i].tor[0]=120;
		}
		if(decstr[i].tor[1]<0) 
			decstr[i].tor[1]=180;
		if(decstr[i].tor[2]<0) 
			decstr[i].tor[2]=290;
		if(decstr[i].len[0]<0) 
			decstr[i].len[0]=float(lencn);
		if(decstr[i].len[1]<0) 
			decstr[i].len[1]=float(lennca);
		if(decstr[i].len[2]<0) 
			decstr[i].len[2]=float(lencac);
		if(decstr[i].ang[0]<0) 
			decstr[i].ang[0]=float(angcacn);
		if(decstr[i].ang[1]<0) 
			decstr[i].ang[1]=float(angcnca);
		if(decstr[i].ang[2]<0) 
			decstr[i].ang[2]=float(angncac);
		//0 1 2 n ca c
		//original phi i-1 psi i-1 omega i
		flagts=tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
			decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
			decstr[i].tor[0]*raddeg,decstr[i].len[0],decstr[i].ang[0]*raddeg,&pn.x,&pn.y,&pn.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinates n %d %f %f %f\n",i,decstr[i].tor[0],decstr[i].len[0],decstr[i].ang[0]);
		}
		decstr[i].ptn.x=pn.x;decstr[i].ptn.y=pn.y;decstr[i].ptn.z=pn.z;
		flagts=tor2pos22(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
			decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
			decstr[i].tor[1]*raddeg,decstr[i].len[1],decstr[i].ang[1]*raddeg,&pt.x,&pt.y,&pt.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinates ca %d %f %f %f\n",i,decstr[i].tor[1],decstr[i].len[1],decstr[i].ang[1]);
		}
		decstr[i].x=pt.x;decstr[i].y=pt.y;decstr[i].z=pt.z;
		flagts=tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
			decstr[i].x,decstr[i].y,decstr[i].z,
			decstr[i].tor[2]*raddeg,decstr[i].len[2],decstr[i].ang[2]*raddeg,&pc.x,&pc.y,&pc.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinates c %d %f %f %f\n",i,decstr[i].tor[2],decstr[i].len[2],decstr[i].ang[2]);
		}
		decstr[i].ptc.x=pc.x;decstr[i].ptc.y=pc.y;decstr[i].ptc.z=pc.z;
		/*
		//atom o
		flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
			decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
			decstr[i].tor[0]*raddeg-PI,lenco,angcaco*raddeg,&pn.x,&pn.y,&pn.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinates o %d %f %f %f\n",i,decstr[i].tor[0],decstr[i].len[0],decstr[i].ang[0]);
		}
		decstr[i-1].pto.x=pn.x;decstr[i-1].pto.y=pn.y;decstr[i-1].pto.z=pn.z;
		*/
		/*
		//atom h
		flagts=bf.tor2pos22(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
			decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,0,lennh,angcnh*raddeg,&pn.x,&pn.y,&pn.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinates d %d\n",i);
		}
		decstr[i].pth.x=pn.x;decstr[i].pth.y=pn.y;decstr[i].pth.z=pn.z;
		*/
	}
	
	/*
	//atom o
	i=seqnum;
	flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
		decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,0,lenco,angcaco*raddeg,&pn.x,&pn.y,&pn.z);
	if(!flagts)
	{
		flagwhole=false;
		printf("wrong front coordinates o %d\n",i);
	}
	decstr[i-1].pto.x=pn.x;decstr[i-1].pto.y=pn.y;decstr[i-1].pto.z=pn.z;
	*/
	/*
	//atom h 
	i=0;
	flagts=bf.tor2pos22(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,decstr[i].x,decstr[i].y,decstr[i].z,
		decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,0,lennh,anghnca*raddeg,&pn.x,&pn.y,&pn.z);
	if(!flagts)
	{
		flagwhole=false;
		printf("wrong front coordinates %d\n",i);
	}
	decstr[i].pth.x=pn.x;decstr[i].pth.y=pn.y;decstr[i].pth.z=pn.z;
	//atom cb
	for(i=0;i<seqnum;i++)
	{
		flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i].x,decstr[i].y,decstr[i].z,
			decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
			torncaccb*raddeg,lenccb,angcaccb*raddeg,&pn.x,&pn.y,&pn.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinates e %d\n",i);
		}
		decstr[i].ptb.x=pn.x;decstr[i].ptb.y=pn.y;decstr[i].ptb.z=pn.z;
		if(decstr[i].aaa=='G')
		{
			decstr[i].ptb.x=decstr[i].x;decstr[i].ptb.y=decstr[i].y;decstr[i].ptb.z=decstr[i].z;
		}
	}
	//atom ha
	for(i=0;i<seqnum;i++)
	{
		flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
			decstr[i].x,decstr[i].y,decstr[i].z,tornccaha*raddeg,lencaha,angccaha*raddeg,&pn.x,&pn.y,&pn.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinates f %d\n",i);
		}
		decstr[i].ptha.x=pn.x;decstr[i].ptha.y=pn.y;decstr[i].ptha.z=pn.z;
	}	
	*/
}
	return flagwhole;
}

bool tor2strp(vector<point3f> &decstr,int seqnum,int istart)
{
	int i;
	point3s pt,pn,pc;
//	BasicFunc bf;
	bool flagts;
	bool flagwhole=true;
	for(i=istart;i<seqnum;i++)
	{
		if(decstr[i].tor[0]<0) 
		{
			if(decstr[i].stype=='E')
				decstr[i].tor[0]=120;
			else if(decstr[i].stype=='H')
				decstr[i].tor[0]=300;
			else
				decstr[i].tor[0]=120;
		}
		if(decstr[i].tor[1]<0) 
			decstr[i].tor[1]=180;
		if(decstr[i].tor[2]<0) 
			decstr[i].tor[2]=290;
		if(decstr[i].len[0]<0) 
			decstr[i].len[0]=float(lencn);
		if(decstr[i].len[1]<0) 
			decstr[i].len[1]=float(lennca);
		if(decstr[i].len[2]<0) 
			decstr[i].len[2]=float(lencac);
		if(decstr[i].ang[0]<0) 
			decstr[i].ang[0]=float(angcacn);
		if(decstr[i].ang[1]<0) 
			decstr[i].ang[1]=float(angcnca);
		if(decstr[i].ang[2]<0) 
			decstr[i].ang[2]=float(angncac);
		//0 1 2 n ca c
		//original phi i-1 psi i-1 omega i
		flagts=tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
			decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
			decstr[i].tor[0]*raddeg,decstr[i].len[0],decstr[i].ang[0]*raddeg,&pn.x,&pn.y,&pn.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinatesp n %d %f %f %f\n",i,decstr[i].tor[0],decstr[i].len[0],decstr[i].ang[0]);
		}
		decstr[i].ptn.x=pn.x;decstr[i].ptn.y=pn.y;decstr[i].ptn.z=pn.z;
		//ca
		flagts=tor2pos22(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
			decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
			decstr[i].tor[1]*raddeg,decstr[i].len[1],decstr[i].ang[1]*raddeg,&pt.x,&pt.y,&pt.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinatesp ca %d %f %f %f\n",i,decstr[i].tor[1],decstr[i].len[1],decstr[i].ang[1]);
		}
		decstr[i].x=pt.x;decstr[i].y=pt.y;decstr[i].z=pt.z;
		//c
		flagts=tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
			decstr[i].x,decstr[i].y,decstr[i].z,
			decstr[i].tor[2]*raddeg,decstr[i].len[2],decstr[i].ang[2]*raddeg,&pc.x,&pc.y,&pc.z);
		if(!flagts)
		{
			flagwhole=false;
			printf("wrong front coordinatesp c %d %f %f %f\n",i,decstr[i].tor[2],decstr[i].len[2],decstr[i].ang[2]);
		}
		decstr[i].ptc.x=pc.x;decstr[i].ptc.y=pc.y;decstr[i].ptc.z=pc.z;
	}
	return flagwhole;
}

int energybondlength(vector<point3f> &decstr,int numseq,double *fene)
{
	int i;
	double totene=0;
//	double meanstd[][4]={
//	{1.466,0.015*4,1.466-0.015*4,1.466+0.015*4},//nca
//	{1.525,0.021*4,1.525-0.021*4,1.525+0.021*4},//cac
//	{1.341,0.016*4,1.341-0.016*4,1.341+0.016*4},//cn
//	{3.813,0.080*4,3.813-0.080*4,3.813+0.080*4},//caca
//	};
	lableproblematic lp1,lp2;
	int istart;
	lp1.nn[3]=0;
	memset(lp1.indn[23],0,numseq*sizeof(int));
	for(i=1;i<numseq;i++)
	{
		if(decstr[i].iaa==5) istart=4;
		else if(decstr[i].iaa==12) istart=7;
		else istart=0;
		if(decstr[i].len[0]<stdlength[istart+2][2])
		{
	//		printf("%3d cn [%.4f %.4f] %.4f\n",i,meanstd[2][2],meanstd[2][3],decstr[i].len[0]);
			lp1.indn[3][lp1.nn[3]]=i;
			lp1.nn[3]++;
			lp1.indn[23][i]=1;
			totene+=stdlength[istart+2][2]-decstr[i].len[0];
		}
		else if(decstr[i].len[0]>stdlength[istart+2][3])
		{
			lp1.indn[3][lp1.nn[3]]=i;
			lp1.nn[3]++;
			lp1.indn[23][i]=1;
			totene+=decstr[i].len[0]-stdlength[istart+2][3];
		}
	}
	lp1.nn[4]=0;
	memset(lp1.indn[24],0,numseq*sizeof(int));
	for(i=0;i<numseq;i++)
	{
		if(decstr[i].iaa==5) istart=4;
		else if(decstr[i].iaa==12) istart=7;
		else istart=0;
		if(decstr[i].len[1]<stdlength[istart+0][2])
		{
	//		printf("%3d nca [%.4f %.4f] %.4f\n",i,meanstd[0][2],meanstd[0][3],decstr[i].len[1]);
			lp1.indn[4][lp1.nn[4]]=i;
			lp1.nn[4]++;
			lp1.indn[24][i]=1;
			totene+=stdlength[istart+0][2]-decstr[i].len[1];
		}
		else if(decstr[i].len[1]>stdlength[istart+0][3])
		{
			lp1.indn[4][lp1.nn[4]]=i;
			lp1.nn[4]++;
			lp1.indn[24][i]=1;
			totene+=decstr[i].len[1]-stdlength[istart+0][3];
		}
	}
	lp1.nn[5]=0;
	memset(lp1.indn[25],0,numseq*sizeof(int));
	for(i=0;i<numseq;i++)
	{
		if(decstr[i].iaa==5) istart=4;
		else if(decstr[i].iaa==12) istart=7;
		else istart=0;
		if(decstr[i].len[2]<stdlength[istart+1][2])
		{
	//		printf("%3d cac [%.4f %.4f] %.4f\n",i,meanstd[1][2],meanstd[1][3],decstr[i].len[2]);
			lp1.indn[5][lp1.nn[5]]=i;
			lp1.nn[5]++;
			lp1.indn[25][i]=1;
			totene+=stdlength[istart+1][2]-decstr[i].len[2];
		}
		else if(decstr[i].len[2]>stdlength[istart+1][3])
		{
			lp1.indn[5][lp1.nn[5]]=i;
			lp1.nn[5]++;
			lp1.indn[25][i]=1;
			totene+=decstr[i].len[2]-stdlength[istart+1][3];
		}
	}
	lp1.nn[6]=0;
	memset(lp1.indn[26],0,numseq*sizeof(int));
	double tdist;
	point3d tp;
	for(i=1;i<numseq;i++)
	{
		tp.x=decstr[i].x-decstr[i-1].x;
		tp.y=decstr[i].y-decstr[i-1].y;
		tp.z=decstr[i].z-decstr[i-1].z;
		tdist=norm(tp);
		if(tdist<stdlength[3][2])
		{
	//		printf("%3d ca [%.4f %.4f] %.4f\n",i,meanstd[3][2],meanstd[3][3],tdist);
			lp1.indn[6][lp1.nn[6]]=i;
			lp1.nn[6]++;
			lp1.indn[26][i]=1;
			totene+=stdlength[3][2]-tdist;
		}
		else if(tdist>stdlength[3][3])
		{
			lp1.indn[6][lp1.nn[6]]=i;
			lp1.nn[6]++;
			lp1.indn[26][i]=1;
			totene+=tdist-stdlength[3][3];
		}
	}
	*fene=totene;
	return lp1.nn[3]+lp1.nn[4]+lp1.nn[5]+lp1.nn[6];
}

string AAtoS(string animoacide)
{
    if(animoacide==" ALA"||animoacide=="AALA"){return ("A") ;}
    if(animoacide==" CYS"||animoacide=="ACYS"){return ("C");}
    if(animoacide==" ASP"||animoacide=="AASP"){return ("D");}
    if(animoacide==" GLU"||animoacide=="AGLU"){return ("E");}
    if(animoacide==" PHE"||animoacide=="APHE"){return ("F");}
    if(animoacide==" GLY"||animoacide=="AGLY"){return ("G");}
    if(animoacide==" HIS"||animoacide=="AHIS"){return ("H");}
    if(animoacide==" ILE"||animoacide=="AILE"){return ("I");}
    if(animoacide==" LYS"||animoacide=="ALYS"){return ("K");}
    if(animoacide==" LEU"||animoacide=="ALEU"){return ("L");}
    if(animoacide==" MET"||animoacide=="AMET"){return ("M");}
    if(animoacide==" ASN"||animoacide=="AASN"){return ("N");}
    if(animoacide==" PRO"||animoacide=="APRO"){return ("P");}
    if(animoacide==" GLN"||animoacide=="AGLN"){return ("Q");}
    if(animoacide==" ARG"||animoacide=="AARG"){return ("R");}
    if(animoacide==" SER"||animoacide=="ASER"){return ("S");}
    if(animoacide==" THR"||animoacide=="ATHR"){return ("T");}
    if(animoacide==" VAL"||animoacide=="AVAL"){return ("V");}
    if(animoacide==" TRP"||animoacide=="ATRP"){return ("W");}
    if(animoacide==" TYR"||animoacide=="ATYR"){return ("Y");}
}

int EMRefinement::fenergybondlength(vector<point3f> &decstr,int numseq,float &fene,bool flagacc)
{
	int i;
	double tdist;
	point3d tp;
//	double meanstd[][4]={
//	{1.466,0.015*4,1.466-0.015*4,1.466+0.015*4},//nca
//	{1.525,0.021*4,1.525-0.021*4,1.525+0.021*4},//cac
//	{1.341,0.016*4,1.341-0.016*4,1.341+0.016*4},//cn
//	{3.813,0.080*4,3.813-0.080*4,3.813+0.080*4},//caca
//	};
	int istart;
if(!flagacc)
{
	lp1.nn[3]=0;
	lp1.fn[3]=0;
	memset(lp1.indn[3],0,numseq*sizeof(int));
	memset(lp1.indf[3],0,numseq*sizeof(double));
	memset(lp1.flap[3],0,numseq*sizeof(bool));
	for(i=1;i<numseq;i++)
	{
		if(decstr[i].resind==5) istart=4;
		else if(decstr[i].resind==12) istart=7;
		else istart=0;
		if(decstr[i].len[0]<stdlength[istart+2][2])
		{
			lp1.nn[3]++;
			lp1.indn[3][i]=1;
			lp1.indf[3][i]=stdlength[istart+2][2]-decstr[i].len[0];
			lp1.fn[3]+=lp1.indf[3][i];
		}
		else if(decstr[i].len[0]>stdlength[istart+2][3])
		{
			lp1.nn[3]++;
			lp1.indn[3][i]=1;
			lp1.indf[3][i]=decstr[i].len[0]-stdlength[istart+2][3];
			lp1.fn[3]+=lp1.indf[3][i];
		}
	}
	lp1.nn[4]=0;
	lp1.fn[4]=0;
	memset(lp1.indn[4],0,numseq*sizeof(int));
	memset(lp1.indf[4],0,numseq*sizeof(double));
	memset(lp1.flap[4],0,numseq*sizeof(bool));
	for(i=0;i<numseq;i++)
	{
		if(decstr[i].resind==5) istart=4;
		else if(decstr[i].resind==12) istart=7;
		else istart=0;
		if(decstr[i].len[1]<stdlength[istart+0][2])
		{
			lp1.nn[4]++;
			lp1.indn[4][i]=1;
			lp1.indf[4][i]=stdlength[istart+0][2]-decstr[i].len[1];
			lp1.fn[4]+=lp1.indf[4][i];
		}
		else if(decstr[i].len[1]>stdlength[istart+0][3])
		{
			lp1.nn[4]++;
			lp1.indn[4][i]=1;
			lp1.indf[4][i]=decstr[i].len[1]-stdlength[istart+0][3];
			lp1.fn[4]+=lp1.indf[4][i];
		}
	}
	lp1.nn[5]=0;
	lp1.fn[5]=0;
	memset(lp1.indn[5],0,numseq*sizeof(int));
	memset(lp1.indf[5],0,numseq*sizeof(double));
	memset(lp1.flap[5],0,numseq*sizeof(bool));
	for(i=0;i<numseq;i++)
	{
		if(decstr[i].resind==5) istart=4;
		else if(decstr[i].resind==12) istart=7;
		else istart=0;
		if(decstr[i].len[2]<stdlength[istart+1][2])
		{
			lp1.nn[5]++;
			lp1.indn[5][i]=1;
			lp1.indf[5][i]=stdlength[istart+1][2]-decstr[i].len[2];
			lp1.fn[5]+=lp1.indf[5][i];
		}
		else if(decstr[i].len[2]>stdlength[istart+1][3])
		{
			lp1.nn[5]++;
			lp1.indn[5][i]=1;
			lp1.indf[5][i]=decstr[i].len[2]-stdlength[istart+1][3];
			lp1.fn[5]+=lp1.indf[5][i];
		}
	}
	lp1.nn[6]=0;
	lp1.fn[6]=0;
	memset(lp1.indn[6],0,numseq*sizeof(int));
	memset(lp1.indf[6],0,numseq*sizeof(double));
	memset(lp1.flap[6],0,numseq*sizeof(bool));
	for(i=1;i<numseq;i++)//i is for d(i-1,i)
	{
		tp.x=decstr[i].ptv[1].x-decstr[i-1].ptv[1].x;
		tp.y=decstr[i].ptv[1].y-decstr[i-1].ptv[1].y;
		tp.z=decstr[i].ptv[1].z-decstr[i-1].ptv[1].z;
		tdist=norm(tp);
		if(tdist<stdlength[3][2])
		{
			lp1.nn[6]++;
			lp1.indn[6][i]=1;
			lp1.indf[6][i]=stdlength[3][2]-tdist;
			lp1.fn[6]+=lp1.indf[6][i];
		}
		else if(tdist>stdlength[3][3])
		{
			lp1.nn[6]++;
			lp1.indn[6][i]=1;
			lp1.indf[6][i]=tdist-stdlength[3][3];
			lp1.fn[6]+=lp1.indf[6][i];
		}
	}
}
else
{
	int iold[4];
	int inew[4];
	double fold[4];
	double fnew[4];
	for(i=0;i<4;i++)
	{
		iold[i]=0;
		inew[i]=0;
		fold[i]=0;
		fnew[i]=0;
	}
	for(i=1;i<numseq;i++) if(lp1.flap[3][i])
	{
		if(decstr[i].resind==5) istart=4;
		else if(decstr[i].resind==12) istart=7;
		else istart=0;
		iold[0]+=lp1.indn[3][i];
		fold[0]+=lp1.indf[3][i];
		if(decstr[i].len[0]<stdlength[istart+2][2])
		{
			lp1.indn[3][i]=1;
			lp1.indf[3][i]=stdlength[istart+2][2]-decstr[i].len[0];
		}
		else if(decstr[i].len[0]>stdlength[istart+2][3])
		{
			lp1.indn[3][i]=1;
			lp1.indf[3][i]=decstr[i].len[0]-stdlength[istart+2][3];
		}
		else
		{
			lp1.indn[3][i]=0;
			lp1.indf[3][i]=0;
		}
		inew[0]+=lp1.indn[3][i];
		fnew[0]+=lp1.indf[3][i];
		lp1.flap[3][i]=false;
	}
	lp1.nn[3]+=inew[0]-iold[0];
	lp1.fn[3]+=fnew[0]-fold[0];

	for(i=0;i<numseq;i++) if(lp1.flap[4][i])
	{
		if(decstr[i].resind==5) istart=4;
		else if(decstr[i].resind==12) istart=7;
		else istart=0;
		iold[1]+=lp1.indn[4][i];
		fold[1]+=lp1.indf[4][i];
		if(decstr[i].len[1]<stdlength[istart+0][2])
		{
			lp1.indn[4][i]=1;
			lp1.indf[4][i]=stdlength[istart+0][2]-decstr[i].len[1];
		}
		else if(decstr[i].len[1]>stdlength[istart+0][3])
		{
			lp1.indn[4][i]=1;
			lp1.indf[4][i]=decstr[i].len[1]-stdlength[istart+0][3];
		}
		else
		{
			lp1.indn[4][i]=0;
			lp1.indf[4][i]=0;
		}
		inew[1]+=lp1.indn[4][i];
		fnew[1]+=lp1.indf[4][i];
		lp1.flap[4][i]=false;
	}
	lp1.nn[4]+=inew[1]-iold[1];
	lp1.fn[4]+=fnew[1]-fold[1];

	for(i=0;i<numseq;i++) if(lp1.flap[5][i])
	{
		if(decstr[i].resind==5) istart=4;
		else if(decstr[i].resind==12) istart=7;
		else istart=0;
		iold[2]+=lp1.indn[5][i];
		fold[2]+=lp1.indf[5][i];
		if(decstr[i].len[2]<stdlength[istart+1][2])
		{
			lp1.indn[5][i]=1;
			lp1.indf[5][i]=stdlength[istart+1][2]-decstr[i].len[2];
		}
		else if(decstr[i].len[2]>stdlength[istart+1][3])
		{
			lp1.indn[5][i]=1;
			lp1.indf[5][i]=decstr[i].len[2]-stdlength[istart+1][3];
		}
		else
		{
			lp1.indn[5][i]=0;
			lp1.indf[5][i]=0;
		}
		inew[2]+=lp1.indn[5][i];
		fnew[2]+=lp1.indf[5][i];
		lp1.flap[5][i]=false;
	}
	lp1.nn[5]+=inew[2]-iold[2];
	lp1.fn[5]+=fnew[2]-fold[2];

	for(i=1;i<numseq;i++) if(lp1.flap[6][i])
	{
		iold[3]+=lp1.indn[6][i];
		fold[3]+=lp1.indf[6][i];
		tp.x=decstr[i].ptv[1].x-decstr[i-1].ptv[1].x;
		tp.y=decstr[i].ptv[1].y-decstr[i-1].ptv[1].y;
		tp.z=decstr[i].ptv[1].z-decstr[i-1].ptv[1].z;
		tdist=norm(tp);
		if(tdist<stdlength[3][2])
		{
			lp1.indn[6][i]=1;
			lp1.indf[6][i]=stdlength[3][2]-tdist;
		}
		else if(tdist>stdlength[3][3])
		{
			lp1.indn[6][i]=1;
			lp1.indf[6][i]=tdist-stdlength[3][3];
		}
		else
		{
			lp1.indn[6][i]=0;
			lp1.indf[6][i]=0;
		}
		inew[3]+=lp1.indn[6][i];
		fnew[3]+=lp1.indf[6][i];
		lp1.flap[6][i]=false;
	}
	lp1.nn[6]+=inew[3]-iold[3];
	lp1.fn[6]+=fnew[3]-fold[3];
}
	*fene=lp1.fn[3]+lp1.fn[4]+lp1.fn[5]+lp1.fn[6];
	return lp1.nn[3]+lp1.nn[4]+lp1.nn[5]+lp1.nn[6];
}