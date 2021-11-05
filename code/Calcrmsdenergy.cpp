double calcrmsdenergy(point3f *decstr,int numseq)
{
    double trmsd=0; 
    double eneval[10];
    int i;
    double enelist[30],enelistbk[30];
    for(i=0;i<20;i++)
    {
        enelist[i]=0;
    }
//    ps.tor2strsg(decstr,numseq,sgpos);
    enelist[0]=energyinitcontact(decstr,numseq);
//  enelist[0]=energypartialrmsd(decstr,numseq);
    trmsd+=enelist[0];
    enelist[1]=0.30*energybondangle(decstr,numseq,&eneval[2]);
    trmsd+=enelist[1]+eneval[2];
    enelist[2]=0.50*energybondlength(decstr,numseq,&eneval[3]);
    trmsd+=enelist[2]+20.0*eneval[3];
//  enelist[3]=energytorsion2(decstr,numseq);
//  trmsd+=enelist[3];
//  calcallenergy(decstr,numseq,enelist);
//  trmsd+=enelist[4];//df
//  trmsd+=enelist[5];//sg
//  trmsd+=enelist[6];//ev
    enelist[3]=energyclash(decstr,numseq,&eneval[4]);
    trmsd+=enelist[3]+eneval[4];
    enelist[4]=4.00*energyrama(decstr,numseq,&eneval[5]);
    trmsd+=enelist[4]+eneval[5];
    enelist[5]=5.00*energyhbondnhoc2(decstr,numseq,&eneval[0],&eneval[1]);
    trmsd+=enelist[5]-1.0*eneval[0]-1.0*eneval[1];
    enelist[6]=0.10*energycontinuoushb(decstr,numseq);
    trmsd+=enelist[6];
    enelist[7]=1.00*energyss2(decstr,numseq);
    trmsd+=enelist[7];
    enelist[8]=1.00*energycontinuoushb2(decstr,numseq);
    trmsd+=enelist[8];
//  printf("ct %6.3f ang %2d len %2d cla %2d hb %2d\n",enelist[0],int(enelist[1]),int(enelist[2]),int(enelist[3]),int(enelist[4]));
    int j;
    for(j=0;j<11;j++) lp1.nn[20+j]=0;
    for(i=0;i<numseq;i++)
    {
        j=lp1.indn[20][i]+lp1.indn[26][i]+lp1.indn[28][i];
        lp1.nn[20]+=j;
        lp1.indf[0][i]=j;
        j=lp1.indn[21][i]+lp1.indn[26][i]+lp1.indn[28][i];
        lp1.nn[21]+=j;
        lp1.indf[1][i]=j;
        if(i<numseq-1) j=lp1.indn[22][i]+lp1.indn[26][i+1]+lp1.indn[28][i];
        else j=lp1.indn[22][i]+lp1.indn[28][i];
        lp1.nn[22]+=j;
        lp1.indf[2][i]=j;
        j=lp1.indn[23][i]+lp1.indn[26][i]+lp1.indn[28][i];
        lp1.nn[23]+=j;
        lp1.indf[3][i]=j;
        j=lp1.indn[24][i]+lp1.indn[26][i]+lp1.indn[28][i];
        lp1.nn[24]+=j;
        lp1.indf[4][i]=j;
        if(i<numseq-1) j=lp1.indn[25][i]+lp1.indn[26][i+1]+lp1.indn[28][i];
        else j=lp1.indn[25][i]+lp1.indn[28][i];
        lp1.nn[25]+=j;
        lp1.indf[5][i]=j;
        j=lp1.indn[26][i]+lp1.indn[28][i];
        lp1.nn[26]+=j;
        lp1.indf[6][i]=j;//ome
        if(i<numseq-1) j=lp1.indn[27][i]+lp1.indn[26][i+1]+lp1.indn[28][i];
        else j=lp1.indn[27][i]+lp1.indn[28][i];
        lp1.nn[27]+=j;
        lp1.indf[7][i]=j;//tor phi psi
        j=lp1.indn[20][i]+lp1.indn[21][i]+lp1.indn[22][i]+lp1.indn[23][i]+lp1.indn[24][i]+lp1.indn[25][i]+
            lp1.indn[26][i]+lp1.indn[27][i]+lp1.indn[28][i];
        lp1.nn[28]+=j;
        lp1.indf[8][i]=j;//sft lmp
        if(i<numseq-1) j=lp1.indn[22][i]+lp1.indn[26][i+1]+lp1.indn[27][i];
        else j=lp1.indn[22][i]+lp1.indn[27][i];
        lp1.nn[29]+=j;
        lp1.indf[9][i]=j;//rot end
        j=lp1.indn[28][i];
        lp1.nn[30]+=j;
        lp1.indf[10][i]=j;//rot middle
    }
    for(j=0;j<11;j++) if(lp1.nn[20+j]!=0)
    {
        for(i=0;i<numseq;i++)
        {
            lp1.indf[j][i]/=double(lp1.nn[20+j]);
        }
        for(i=1;i<numseq;i++)
        {
            lp1.indf[j][i]+=lp1.indf[j][i-1];
        }
        for(i=numseq-1;i>0;i--)
        {
            lp1.indf[j][i]=lp1.indf[j][i-1];
        }
        lp1.indf[j][0]=0;
    }
    return trmsd;
}

double energyinitcontact(point3f *decstr,int numseq)
{
    if(npaa3==0) return 0.0;
    double totene=0;
    int i,j,k;
    point3d tp;
    double tval;
//    pairaa* paa3;
    for(i=0;i<numseq;i++)
    {
        for(j=i+2;j<numseq;j++)
        {
            k=j*numseq+i;
            tp.x=decstr[j].x-decstr[i].x;
            tp.y=decstr[j].y-decstr[i].y;
            tp.z=decstr[j].z-decstr[i].z;
            paa3[k].pval=bf.norm(tp);
        }
    }
    for(i=0;i<npaa3;i++)
    {
        k=paa3[i].tnum*numseq+paa3[i].pbin;
    //  tval=(paa3[k].pval-paa3[i].dist)*(paa3[k].pval-paa3[i].dist);//deltad^2;
    //  tval=(paa3[k].pval-paa3[i].dist)*(paa3[k].pval-paa3[i].dist)/paa3[i].dstd/paa3[i].dstd;//deltad^2
    //  tval=fabs(paa3[k].pval-paa3[i].dist);//deltad;
        tval=fabs(paa3[k].pval-paa3[i].dist)/paa3[i].dstd;//deltad
    //  tval=-exp(-(paa3[k].pval-paa3[i].dist)*(paa3[k].pval-paa3[i].dist));//-exp(-deltad^2)
    //  tval=-exp(-(paa3[k].pval-paa3[i].dist)*(paa3[k].pval-paa3[i].dist)/paa3[i].dstd/paa3[i].dstd);//-exp(-deltad^2)
        totene+=tval;
    }
    return contwt*totene/double(npaa3);
}
