#include "ppin.h"

void PrintCorrelation(xsimilar,size)
float    **xsimilar;
int        size;
{
 int        i,ii,j,jj,q;
 float      maxsim,minsim,meansim,maxsim2,minsim2,meansim2,rmsdsim,rmsdsim2;
 float      h,maxnn,step2,q11,q12,q13,q14,q15,q16,q21,q22,q23,q24,q25,q26;
 
         printf("STOP forced to check Matrix of Similarity values\n");
         maxsim=minsim=xsimilar[0][1];
         maxsim2=minsim2=fabs(xsimilar[0][1]);
         meansim=meansim2=0.0;
         rmsdsim=rmsdsim2=0.0;
         ii=0;
         for (i=0;i<size;i++){ for (j=i+1;j<size;j++){
             if (xsimilar[i][j]>maxsim)maxsim=xsimilar[i][j];
             if (xsimilar[i][j]<minsim)minsim=xsimilar[i][j];
             if (fabs(xsimilar[i][j])>maxsim2)maxsim2=fabs(xsimilar[i][j]);
             if (fabs(xsimilar[i][j])<minsim2)minsim2=fabs(xsimilar[i][j]);
             meansim+=xsimilar[i][j];
             meansim2+=fabs(xsimilar[i][j]);
             rmsdsim+=xsimilar[i][j]*xsimilar[i][j];
             rmsdsim2+=xsimilar[i][j]*xsimilar[i][j];
             ii++;
             }}
         if (ii>0)meansim=meansim/ii;
         if (ii>0)meansim2=meansim2/ii;
         if (ii>0)rmsdsim=rmsdsim/ii;
         if (ii>0)rmsdsim2=rmsdsim2/ii;
         rmsdsim2=rmsdsim2-meansim2*meansim2;
         rmsdsim=rmsdsim-meansim*meansim;
         if (rmsdsim>0)rmsdsim=sqrt(rmsdsim);
         if (rmsdsim2>0)rmsdsim2=sqrt(rmsdsim2);
         if (ii>0) maxnn=1+(log(ii)/log(2));
         else      maxnn=100;
         h=(maxsim-minsim)/maxnn;
         q11=q12=q13=q14=q15=q16=0.0;
         q=0;
         printf ("\nQuartil analysis of Distance Similarity between Genes (k=%d h=%f)\n",(int)maxnn,h); 
         for (jj=0;jj<maxnn;jj++){
           step2=minsim+h*(jj+1);
           printf("Threshold %10.5f \t#Edges %d\t Percentage %10.2f \n",step2,q,100*(float)q/ii);
           q=0;
           for (i=0;i<size;i++) for (j=i+1;j<size;j++) if (xsimilar[i][j]<step2) q++;
           if ((float) q/ii > 0.25 && q11==0) q11=step2;  
           if ((float) q/ii > 0.50 && q12==0) q12=step2;  
           if ((float) q/ii > 0.75 && q13==0) q13=step2;
           if ((float) q/ii > 0.90 && q14==0) q14=step2;
           if ((float) q/ii > 0.95 && q15==0) q15=step2;
           if ((float) q/ii > 0.99 && q16==0) q16=step2;
         }
         
         h=(maxsim2-minsim2)/maxnn;
         printf ("\nQuartil analysis of Absolute Distance Similarity between Genes (k=%d h=%f)\n",(int)maxnn,h); 
         q21=q22=q23=q24=q25=q26=0.0;
         q=0;
         for (jj=0;jj<maxnn;jj++){
           step2=minsim2+h*(jj+1);
           //printf("step2=%f q21=%f q22=%f q23=%f q24=%f q25=%f q26=%f Q=%d II=%d Q/II=%f \n",step2,q21,q22,q23,q24,q25,q26,q,ii,(float)q/ii);
           printf("Threshold %10.5f \t#Edges %d\t Percentage %10.2f \n",step2,q,100*(float)q/ii);
           q=0;
           for (i=0;i<size;i++) for (j=i+1;j<size;j++) if (fabs(xsimilar[i][j])<step2) q++;
           if ((float) q/ii > 0.25 && q21==0) q21=step2;  
           if ((float) q/ii > 0.50 && q22==0) q22=step2;  
           if ((float) q/ii > 0.75 && q23==0) q23=step2;
           if ((float) q/ii > 0.90 && q24==0) q24=step2;
           if ((float) q/ii > 0.95 && q25==0) q25=step2;
           if ((float) q/ii > 0.99 && q26==0) q26=step2;
         }
   

         printf("\nMaximum Distance Similarity between Genes: %f\n",maxsim);
         printf("Minimum Distance Similarity between Genes: %f\n",minsim);
         printf("Average Distance Similarity between Genes: %f\n",meansim);
         printf("RMSD    Distance Similarity between Genes: %f\n",rmsdsim);
         printf("Quartil Distance Similarity between Genes: %f (Q1) \t %f (Q2) \t %f (Q3)\n", q11,q12,q13);
         printf("Top 10  Distance Similarity between Genes: %f (90%)\t%f (95%)\t%f (99%) \n", q14,q15,q16);
         printf("Maximum Absolute Distance Similarity between Genes: %f\n",maxsim2);
         printf("Minimum Absolute Distance Similarity between Genes: %f\n",minsim2);
         printf("Average Absolute Distance Similarity between Genes: %f\n",meansim2);
         printf("RMSD    Absolute Distance Similarity between Genes: %f\n",rmsdsim2);
         printf("Quartil Absolute Distance Similarity between Genes: %f (Q1) \t %f (Q2) \t %f (Q3) \n", q21,q22,q23);
         printf("Top 10  Absolute Distance Similarity between Genes: %f (90%)\t%f (95%)\t%f (99%) \n", q24,q25,q26);
         exit(1);
}
