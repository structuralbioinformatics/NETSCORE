#include "ppin.h"

void XRankProd(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples)
float    **xPvalue,**xES;
express  *xprobe;
express  *xZscore;
int      **xgroup,*xdgroup,xngroup,xPhenoRef,*xmethod,number_phe,number_samples[];

{
   int   i,j,k,n,m,s,ii,iii,jj,kk,nr,nn,n1,n2,mm,ss,sss,mRNA,genes,n_probe,irank,*skip,*gene_probe,*fragment_probe,*ig_up,*ig_down;
   int   dimension_samples,dimension_probes,log_scale,count,count_ref,ii_probe,ii_found,iup,idown,x,y,z,top;
   int   *ivector();
   float a,b,emax,emin,fold,pmax,mean,std,rmsd,median,r_up,r_down;
   float **matrix(),*fvector(),*generate_random();
   float **data,*sorted,*gene,*g_up,*g_down,*random_rp, *fc,*e_up, *e_down, *f_down, *f_up,  *rp_down,*rp_up;
   void  sort2();

   nr=100; //number of random samples
   fold=(float)(xmethod[19])/1000;
   top = xmethod[20];
   pmax=(float)(xmethod[21]/1000000); 
  

   if (fold<0) {fold=1.5;}
   if (pmax<=0) {pmax=1.0e-6;}
   nn=mm=0;
   nn=number_samples[xPhenoRef];
   for (j=0;j<number_phe;j++){
       if (j!=xPhenoRef){mm+=number_samples[j];}}
   dimension_samples=nn*mm;  
  
   dimension_probes=0;
   for (i=0;i<xngroup;i++){
   for (n=0;n<xdgroup[i];n++){   
       dimension_probes+=xprobe[xgroup[i][n]].split;
   }}
   for (i=0;i<xngroup;i++){
   for (n=0;n<xdgroup[i];n++){   
       ii=xgroup[i][n];
       if (ii>dimension_probes)dimension_probes=ii;
   }}

   log_scale=0;
   for (i=0;i<xngroup;i++){
   for (n=0;n<xdgroup[i];n++){
   for (j=0;j<number_phe;j++){
   for (s=0;s<number_samples[j];s++){
     if (xprobe[xgroup[i][n]].fragment[kk].phenotype[j].defined[s]>0.0){
     if (xprobe[xgroup[i][n]].fragment[kk].phenotype[xPhenoRef].sample[s]<0){
        log_scale=1;
     }}}}}}

   if (log_scale==1){printf("Execute XRankProd with logarithm scale\n");}
   else             {printf("Execute XRankProd with raw expression\n");}

   data=matrix(0,dimension_samples+1,0,dimension_probes+1);
   skip=ivector(0,dimension_probes+1);
   gene_probe=ivector(0,dimension_probes+1);
   fragment_probe=ivector(0,dimension_probes+1);
   

   for (i=0;i<dimension_probes+1;i++){
      skip[i]=0;
      gene_probe[i]=0;
      for (j=0;j<dimension_samples+1;j++){
        if (log_scale==1){
          data[j][i]=0.0;
        }else{
          data[j][i]=1.0;
        }
       }}

   mRNA=0;
   genes=0;
   n_probe=0;
   emax=0.0; 
   emin=0.0;
  
   if (xmethod[4]==0){
    for (i=0;i<xngroup;i++){
    for (n=0;n<xdgroup[i];n++){
     ii=xgroup[i][n];
     if (ii>genes){genes=ii;}
     if (skip[ii]==0){
       skip[ii]=1;
       for (s=0;s<number_samples[xPhenoRef];s++){
         n1=0;
         a=0.0;
         for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
           if (xprobe[xgroup[i][n]].fragment[kk].phenotype[xPhenoRef].defined[s]>0.0){
              a+=xprobe[xgroup[i][n]].fragment[kk].phenotype[xPhenoRef].sample[s];
              n1++;
           }}
         if (n1>0){ a=a/n1;}
         if (a>emax){emax=a;}
         if (emin==0.0 || a<emin){emin=a;}
         sss=0;
         for (j=0;j<number_phe;j++){
           if (j!=xPhenoRef){
             for (ss=0;ss<number_samples[j];ss++){
               n2=0;
               b=0.0;
               for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
                 if (xprobe[xgroup[i][n]].fragment[kk].phenotype[j].defined[ss]>0.0){
                    b+=xprobe[xgroup[i][n]].fragment[kk].phenotype[j].sample[ss];
                    n2++;
                 }}  
               if (n2>0) {b=b/n2;}
               if (b>emax){emax=b;}
               if (emin==0.0 || b<emin){emin=b;}
               jj=mm*s + sss ;
               sss++;
               if (n1>0 && n2>0 && b!=0.0 && a!=0.0) {
                 if (log_scale==0){data[jj][ii]=a/b;}else{data[jj][ii]=a-b;}
                 //printf("DATA[%d][%d]= %f f(%f,%f) NN %d MM %d S %d SS %d SSS %d (I=%d N=%d) MAXSIZE %d %d\n",jj,ii,data[jj][ii],a,b,nn,mm,s,ss,sss,i,n,dimension_samples,dimension_probes);
                 }
              }
           }}
        }
      }}}
      genes = genes +1;
      mRNA  = genes;
   }else{
    for (i=0;i<xngroup;i++){
    for (n=0;n<xdgroup[i];n++){
      ii=xgroup[i][n];
      if (skip[ii]==0){
        skip[ii]=1;
        a=b=0.0;
        for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
         for (s=0;s<number_samples[xPhenoRef];s++){
             if (xprobe[ii].fragment[kk].phenotype[xPhenoRef].defined[s]>0.0){ a=xprobe[ii].fragment[kk].phenotype[xPhenoRef].sample[s]; }
             sss=0;
             if (a>emax){emax=a;}
             if (emin==0.0 || a<emin){emin=a;}
             for (j=0;j<number_phe;j++){
               if (j!=xPhenoRef){
                 for (ss=0;ss<number_samples[j];ss++){
                   if (xprobe[ii].fragment[kk].phenotype[j].defined[ss]>0.0){ b=xprobe[ii].fragment[kk].phenotype[j].sample[ss];}
                   jj=mm*s + sss ;
                   sss++;
                   if (b>emax){emax=b;}
                   if (emin==0.0 || b<emin){emin=b;}
                   if (b!=0.0 && a!=0.0) {if (log_scale==0){data[jj][n_probe]=a/b;}else{data[jj][n_probe]=a-b;}}
                   //printf("DATA[%d][%d]= %f f(%f,%f) NN %d MM %d S %d SS %d SSS %d (I=%d N=%d) MAXSIZE %d %d\n",jj,n_probe,data[jj][n_probe],a,b,nn,mm,s,ss,sss,i,n,dimension_samples,dimension_probes);
                }}}
           }
         gene_probe[n_probe]=ii;
         fragment_probe[n_probe]=kk;
         n_probe++;
         if (n_probe>dimension_probes){printf("ERROR[%d,%d,%d] n_probes %d dimension %d\n",i,n,ii,n_probe,dimension_probes);}
        }
        if (xprobe[xgroup[i][n]].split==0){
         gene_probe[n_probe]=ii;
         fragment_probe[n_probe]=0;
         n_probe++;
         if (n_probe>dimension_probes){printf("ERROR[%d,%d,%d] n_probes %d dimension %d\n",i,n,ii,n_probe,dimension_probes);}
        }
     }}}
     mRNA=n_probe;
   }

   if (top==0){top=mRNA;}

    sorted =fvector(0,mRNA+1);
    rp_up  =fvector(0,mRNA+1);
    rp_down=fvector(0,mRNA+1);
    fc     =fvector(0,mRNA+1);
    e_up   =fvector(0,mRNA+1);
    e_down =fvector(0,mRNA+1);
    f_up   =fvector(0,mRNA+1);
    f_down =fvector(0,mRNA+1);
    gene   =fvector(0,mRNA+1);
    g_up   =fvector(0,mRNA+1);
    g_down =fvector(0,mRNA+1);
    ig_up  =ivector(0,mRNA+1);
    ig_down=ivector(0,mRNA+1);
    random_rp=fvector(0,mRNA*nr+1);

    for (i=0;i<mRNA+1;i++){
        rp_up[i]=0.0;
        rp_down[i]=0.0;
        fc[i]=0.0;
        e_up[i]=mRNA;
        e_down[i]=mRNA;
        f_up[i]=1.0;
        f_down[i]=1.0;
        gene[i]=(float)(i);
       }
    printf("XRankProd ranking mRNA\n");
    for (s=0;s<dimension_samples;s++){
       for (i=1;i<mRNA+1;i++){ sorted[i]=data[s][i-1];gene[i]=(float)(i);}
       //if (s==0){ for (i=1;i<mRNA+1;i++){printf("S %d SORTED[%d]=%f GENE[%d]= %f\n",s,i,sorted[i],i,gene[i]);}}
       sort2(mRNA,sorted,gene);
       //if (s==0){ for (i=1;i<mRNA+1;i++){printf("S %d NEW SORTED[%d]=%f GENE[%d]= %f\n",s,i,sorted[i],i,gene[i]);}}
       for (j=mRNA;j>0;j--){
         i=mRNA -j + 1;
         ii=(int)(gene[j]);
         rp_up[ii]    += log10f((float)(i)/mRNA);
         g_up[ii]      = (float)(ii); 
         rp_down[ii]  += log10f((float)(j)/mRNA);  
         g_down[ii]    = (float)(ii); 
         fc[ii]       += sorted[j]/dimension_samples;
         }
       }
/*
    for (j=mRNA;j>0;j--){
      ii=(int)(gene[j]);
      printf("REAL Probe %d  RP_UP[%d]=%f RP_DOWN[%d]= %f\n",j,ii,rp_up[ii],ii,rp_down[ii]);
      }
*/
    mean=std=0.0;
    for (s=1;s<mRNA+1;s++){
       mean+=rp_up[s]/mRNA;
       std +=rp_up[s]*rp_up[s]/mRNA;
    }
    if (std>mean*mean){rmsd=sqrtf(std-mean*mean);}else{rmsd=.0;}
    printf("Distribution of RP_UP  (Mean %f rmsd %f Q %f Median %f)\n",mean,rmsd,std,rp_up[(int)(mRNA/2)]);
    mean=std=0.0;
    for (s=1;s<mRNA+1;s++){
       mean+=rp_down[s]/mRNA;
       std +=rp_down[s]*rp_down[s]/mRNA;
    }
    if (std>mean*mean){rmsd=sqrtf(std-mean*mean);}else{rmsd=.0;}
    printf("Distribution of RP_DOWN  (Mean %f rmsd %f Q %f Median %f)\n",mean,rmsd,std,rp_down[(int)(mRNA/2)]);

    printf("Generate Random XRankProd\n");
    random_rp=generate_random(nr,mRNA,number_phe,number_samples,xPhenoRef,emax,emin,log_scale);
    mean=std=0.0;
    for (s=1;s<mRNA*nr+1;s++){
       mean+=random_rp[s]/(mRNA*nr);
       std +=random_rp[s]*random_rp[s]/(mRNA*nr);
    }
    if (std>mean*mean){rmsd=sqrtf(std-mean*mean);}else{rmsd=.0;}
    printf("Calculate p-values for XRankProd (Mean %f rmsd %f Q %f Median %f)\n",mean,rmsd,std,random_rp[(int)(mRNA*nr/2)]);
    
    sort2(mRNA,rp_up,g_up);
    for (j=mRNA;j>0;j--){
       ii=(int)(g_up[j]);
       ig_up[ii]=j;
       if ( ig_up[ii] < top) {printf("REORDER UP %d < %d Probe %d RP_UP[%d]= %f RP_DOWN[%d]= %f FC[%d]= %f\n",j,top,ii ,gene_probe[ii-1],rp_up[j],gene_probe[ii-1],rp_down[j],ii,fc[ii]);}
    }
    x=1;
    y=0;
    z=0;
    for (i=mRNA;i>0;i--){e_up[i]=f_up[i]=0.0;}
    for (i=1;i<mRNA+1;i++){
    //for (i=mRNA;i>0;i--){
      y++;
      ii=(int)(g_up[i]);
      for (s=x;s<mRNA*nr+1;s++){
      //for (s=mRNA*nr;s>x;s--){
         if (random_rp[s]==0){z=s;}
         if (random_rp[s]>rp_up[i] && random_rp[s]!=0){
          x=s;
          e_up[i]=(float)(s-z-1)/nr;
          f_up[i]=e_up[i]/y;
          //printf("UP S= %d mRNA*NR= %d R_UP[%d in %d = %d]= %f Random[%d]= %f X= %d Y= %d Probe %d FC[%d]= %f Freq_up[%d]= %f E.up[%d]= %f\n",s,mRNA*nr,ii,ig_up[ii],i,rp_up[ig_up[ii]],s,random_rp[s], x, y, (int)(g_up[i]), gene_probe[ii-1], fc[ii],i,f_up[i],i,e_up[i]);
          break;
          }
      }}

    sort2(mRNA,rp_down,g_down);
    for (j=mRNA;j>0;j--){
       ii=(int)(g_down[j]);
       ig_down[ii]=j;
       if (ig_down[ii]<top){printf("REORDER DOWN %d < %d Probe %d RP_UP[%d]= %f RP_DOWN[%d]= %f FC[%d]= %f\n",j,top,ii,gene_probe[ii],rp_up[j],gene_probe[ii],rp_down[j],ii,fc[ii]);}
    }
    x=1;
    y=0;
    z=0;
    for (i=mRNA;i>0;i--){e_down[i]=f_down[i]=0.0;}
    for (i=1;i<mRNA+1;i++){
    //for (i=mRNA;i>0;i--){
      y++;
      for (s=x;s<mRNA*nr+1;s++){
      //for (s=mRNA*nr;s>x;s--){
         if (random_rp[s]==0){z=s;}
         if (random_rp[s]>rp_down[i] && random_rp[s]!=0){
          x=s;
          e_down[i]=(float)(s-z-1)/nr;
          f_down[i]=e_down[i]/y;
          //printf("DOWN S= %d mRNA*NR= %d R_DOWN[%d in %d = %d]= %f Random[%d]= %f X= %d Y= %d Probe %d FC[%d]= %f Freq_up[%d]= %f E.up[%d]= %f\n",s,mRNA*nr,ii,ig_up[ii],i,rp_down[i],s,random_rp[s], x, y, (int)(g_up[i]), gene_probe[ii-1], fc[ii],i,f_up[i],i,e_up[i]);
          break;
          }
      }}

  
   if (xmethod[4]==0){
     for (i=0;i<xngroup;i++){
     for (n=0;n<xdgroup[i];n++){
      ii=xgroup[i][n];
      iii=ii+1;
      if (iii<=mRNA){
       iup  =ig_up[iii];
       idown=ig_down[iii];
       printf("Calculate xES[%d][%d] UP %d (%e) DOWN %d (%e) FC %e\n",ii,xPhenoRef,iup,idown,f_up[iup],f_down[idown],fc[iii]);
       count=count_ref=0;
       for (j=0;j<number_phe;j++){
       for (s=0;s<number_samples[j];s++){
       for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
        if (xprobe[xgroup[i][n]].fragment[kk].phenotype[j].defined[s]>0.0){ 
          if (j==xPhenoRef){count_ref++;}else{count++;}
         }
        }}}
       if (count>0 && count_ref>0 && (iup<top||idown<top)&&(f_up[iup]<pmax||f_down[idown]<pmax)){
        if (fc[iii] > fold && (f_up[iup]<pmax||f_down[idown]<pmax) ){
         xES[ii][xPhenoRef]=1.0;
        }else{
         if (log_scale==0 && fc[iii] < 1/fold){xES[ii][xPhenoRef]=1.0;}
         if (log_scale==1 && fc[iii] < -fold){xES[ii][xPhenoRef]=1.0;}
        }    
        if ( (log_scale==0 && fc[iii]>1 ) ||  (log_scale==1 && fc[iii]>0 ) ){
         xPvalue[ii][xPhenoRef] = f_up[iup];
        }else{
         xPvalue[ii][xPhenoRef] = f_down[idown];
        }
       }
       if (xES[ii][xPhenoRef]>0){printf(" -- xES[%d][%d]= %f FC=%f RP-up=%f RP-down=%f  P-value=%e E-Value-Down= %f E-Value-Up= %f\n",ii,xPhenoRef,xES[ii][xPhenoRef],fc[iii], rp_up[iup],rp_down[idown],xPvalue[ii][xPhenoRef],e_down[idown],e_up[iup]);}
     }}}
   }else{
     for (iii=1;iii<mRNA+1;iii++){
      n_probe = iii-1;
      ii   = gene_probe[n_probe];
      iup  =ig_up[iii];
      idown=ig_down[iii];
      count=count_ref=0;
      if (iii==1 || ii!=ii_found){r_up=r_down=0.0;}
      ii_probe=fragment_probe[n_probe];
      for (j=0;j<number_phe;j++){
      for (s=0;s<number_samples[j];s++){
      for (kk=0;kk<xprobe[ii].split;kk++){
       if (xprobe[ii].fragment[kk].phenotype[j].defined[s]>0.0){ 
         if (j==xPhenoRef){count_ref++;}else{count++;}
        }
       }}}
      if (count>0 && count_ref>0 && (iup<top||idown<top)&&(f_up[iup]<pmax||f_down[idown]<pmax)){
      if (
          xES[ii][xPhenoRef]==0.0 || 
         (xES[ii][xPhenoRef]==1.0 && ((log_scale==0 && fc[iii]>1 ) || (log_scale==1 && fc[iii]>0 ))  && (xPvalue[ii][xPhenoRef]>f_up[iup] || r_up>rp_up[iup]) ) || 
         (xES[ii][xPhenoRef]==1.0 && ((log_scale==0 && fc[iii]<1 ) || (log_scale==1 && fc[iii]<0 ))  && (xPvalue[ii][xPhenoRef]>f_down[idown] || r_down> rp_down[idown]) ) 
         ){ 
       if (fc[iii] > fold ){
        xES[ii][xPhenoRef]=1.0;
       }else{
        if (log_scale==0 && fc[iii] < 1/fold){xES[ii][xPhenoRef]=1.0;}
        if (log_scale==1 && fc[iii] < -fold){xES[ii][xPhenoRef]=1.0;}
       }   
       if  ( xES[ii][xPhenoRef]==1.0){
        ii_found=ii;
        r_up= rp_up[iup];
        r_down= rp_down[idown];
        if ( (log_scale==0 && fc[iii]>1 ) ||  (log_scale==1 && fc[iii]>0 ) ){
         xPvalue[ii][xPhenoRef] = f_up[iup];
        }else{
         xPvalue[ii][xPhenoRef] = f_down[idown];
        }}
        if (xES[ii][xPhenoRef]>0){printf("Selected  Probe[%d = %d] %s xES[%d][%d]= %f FC=%f RP-up=%f (use %f) RP-down=%f (use %f) P-value=%e E-Value-Down= %f E-Value-Up= %f f_up[%d]= %f  f_down[%d]= %f\n",iii,ii_probe,xprobe[ii].fragment[ii_probe].name,ii,xPhenoRef,xES[ii][xPhenoRef],fc[iii], rp_up[iup],r_up,rp_down[idown],r_down, xPvalue[ii][xPhenoRef],e_down[idown],e_up[iup],iup,f_up[iup],idown,f_down[idown]);}
      }else{
        if (xES[ii][xPhenoRef]>0){printf("Discarded Probe[%d = %d] %s xES[%d][%d]= %f FC=%f RP-up=%f (use %f) RP-down=%f (use %f) P-value=%e E-Value-Down= %f E-Value-Up= %f f_up[%d]= %f  f_down[%d]= %f\n",iii,ii_probe,xprobe[ii].fragment[ii_probe].name,ii,xPhenoRef,xES[ii][xPhenoRef],fc[iii], rp_up[iup],r_up,rp_down[idown],r_down, xPvalue[ii][xPhenoRef],e_down[idown],e_up[iup],iup,f_up[iup],idown,f_down[idown]);}
      }}
     }
   }  
 
   printf("XRankProd Clean up memory ( genes %d probes %d samples %d genes*nr %d) \n",mRNA,dimension_probes,dimension_samples,mRNA*nr);
   printf("XRankProd Clean up SORT %d\n");
   free_fvector(sorted,0,mRNA+1);
   printf("XRankProd Clean up RP\n");
   free_fvector(rp_up,0,mRNA+1);
   free_fvector(rp_down,0,mRNA+1);
   printf("XRankProd Clean up FC\n");
   free_fvector(fc,0,mRNA+1);
   printf("XRankProd Clean up E\n");
   free_fvector(e_up,0,mRNA+1);
   free_fvector(e_down,0,mRNA+1);
   printf("XRankProd Clean up F\n");
   free_fvector(f_up,0,mRNA+1);
   free_fvector(f_down,0,mRNA+1);
   printf("XRankProd Clean up GENE\n");
   free_fvector(gene,0,mRNA+1);
   free_fvector(g_up,0,mRNA+1);
   free_fvector(g_down,0,mRNA+1);
   free_ivector(ig_up,0,mRNA+1);
   free_ivector(ig_down,0,mRNA+1);
   printf("XRankProd Clean up SKIP %d\n",dimension_probes);
   free_ivector(skip,0,dimension_probes+1);
   printf("XRankProd Clean up DATA [%d]x[%d]\n",dimension_samples,dimension_probes);
   free_matrix(data,0,dimension_samples+1,0,dimension_probes+1);
   printf("XRankProd Clean up RANDOM %d\n",mRNA*nr);
   free_fvector(random_rp,0,mRNA*nr+1);
   printf("XRankProd Clean up PROBES %d\n",dimension_probes);
   free_ivector(gene_probe,0,dimension_probes+1);
   free_ivector(fragment_probe,0,dimension_probes+1);
   printf("Done XRankProd \n");
       
}    
        
float *generate_random(nr,mRNA,number_phe,number_samples,xPhenoRef,emax,emin,log_scale)
int nr,mRNA,xPhenoRef,number_phe,*number_samples,log_scale;
float  emax,emin;
{
 int   i,j,jj,ii,ss,k,kk,s,sss,dim,iran,trial,nn,mm;
 float *rp, *all_rp, rand_a, rand_b, **rand_data,*sorted,*gene,value;
 float *fvector(),**matrix(),ran3();
 void sort(),sort2();

 nn=mm=0;
 nn=number_samples[xPhenoRef];
 for (j=0;j<number_phe;j++){
       if (j!=xPhenoRef){mm+=number_samples[j];}}
 dim=nn*mm;  

 sorted =fvector(0,mRNA+1);
 gene   =fvector(0,mRNA+1);
 rand_data=matrix(0,dim+1,0,mRNA+1);
 rp     =fvector(0,mRNA+1);
 all_rp =fvector(0,mRNA*nr+1);

 for (j=0;j<mRNA+1;j++){ rp[j]=0.0;}
 for (kk=0;kk<nr;kk++){
 for (j=0;j<mRNA+1;j++){
    all_rp[kk*mRNA+j]=rp[j];
   }}
   
 for (k=0;k<nr;k++){

   for (j=0;j<mRNA+1;j++){ rp[j]=0.0;}

   for (i=0;i<mRNA+1;i++){
   for (j=0;j<dim;j++){
      if (log_scale==1){
          rand_data[j][i]=0.0;
      }else{
          rand_data[j][i]=1.0;
      }
   }}

   for (ii=1;ii<mRNA+1;ii++){ 
    gene[ii]=(float)(ii);
    for (s=0;s<number_samples[xPhenoRef];s++){
         iran=(int)(ran3(&s));
         trial=0;
         if (log_scale==0){value=emax*fabs(ran3(&iran));}else{value=emax*ran3(&iran);}
         while (trial<100){
            trial++;
            if (value>=emin){ rand_a=value;break;}
            }
         sss=0;
         for (j=0;j<number_phe;j++){
            if (j!=xPhenoRef){
            for (ss=0;ss<number_samples[j];ss++){
                   jj=mm*s + sss ;
                   sss++;
                   iran=(int)(ran3(&ss));
                   trial=0;
                   if (log_scale==0){value=emax*fabs(ran3(&iran));}else{value=emax*ran3(&iran);}
                   while (trial<100){
                     trial++;
                     if (value>=emin){ rand_b=value;break;}
                     }
                   if (log_scale==0){
                     if (rand_a>0 &&  rand_b>0){rand_data[jj][ii]= rand_a / rand_b;}
                   }else{
                    rand_data[jj][ii]= rand_a - rand_b;
                   }
                   //printf("RAN_DATA[%d][%d]= %f f(%f,%f) NN %d MM %d S %d SS %d SSS %d \n",jj,ii,rand_data[jj][ii],rand_a,rand_b,nn,mm,s,ss,sss);
            }}
            }
         }
    } 
   for (s=0;s<dim;s++){
       for (i=1;i<mRNA+1;i++){ sorted[i]=rand_data[s][i];gene[i]=(float)(i);}
       sort2(mRNA,sorted,gene);
       for (j=mRNA;j>0;j--){
         i=mRNA -j + 1;
         ii=(int)(gene[j]);
         rp[ii]    += log10f((float)(i)/mRNA);
         }
       }
   
/*
   for (j=mRNA;j>0;j--){
      ii=(int)(gene[j]);
      printf("RANDOM Probe %d  RP[%d]=%f \n",j,ii,rp[ii]);
      }
*/

   for (j=mRNA;j>0;j--){
    all_rp[k*mRNA+j]=rp[j];
   }
 }
 sort(mRNA*nr,all_rp);
 

 free_fvector(sorted,0,mRNA+1);
 free_fvector(gene,0,mRNA+1);
 free_fvector(rp,0,mRNA+1);
 free_matrix(rand_data,0,dim+1,0,mRNA+1);

 return all_rp;

}


