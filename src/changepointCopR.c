#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>

/****************************************************************************************/
/*   This file contains the functions to calculate the change-point test for copulas    */
/*                                                                                      */
/* par: Bruno Remillard, le 11 mai  2010                                                */
/*                                                                                      */
/****************************************************************************************/

   double maxi(double u, double v)

   {
      if( u> v)
         return u;
      else
         return v;

   }

   double mini(double u, double v)

   {
      if( u> v)
         return v;
      else
         return u;

   }





   void rank(double *x, double *r, int n)

   {



      int i, j;
      int count;

      for(i=0;i<n;i++)
      {
         count=0;
         for(j=0;j<n;j++)
         {
            if(x[j] <= x[i])
               ++count;
         }
         r[i] = (double)count;

      }
   }



   double mean(double *x, int n)

   {
      int i;
      double sum = 0.0;

      for(i=0;i<n;i++)
         sum += x[i];

      return sum/((double) n);
   }



   double sum(double *x, int n)

   {
      int i;
      double s = 0.0;

      for(i=0;i<n;i++)
         s += x[i];

      return s;
   }


   double maxvec(double *x, int n)

   {
      int i;
      double y, s;

      s = 0.0;
      for(i=0;i<n;i++)
      {
         y = fabs(x[i]);

         if(s <y)
            s = y;
        /* printf("s = %f\n",s); */


      }

      return s;
   }


   void multvec(double *x, double *y, double *xy, int n)

   {
      int i;


      for(i=0;i<n;i++)
         xy[i] = x[i]*y[i];


   }





   void empcdf(double *x, int *n, int *d, double *u, double *M, double *cumsum)
   {
       int i,j;
       int prod, sum = 0;
       double temp, p;




      // printf("n = %d d = %d\n", n[0],d[0]);



       for(i=0;i<n[0];i++)
       {
           prod = 1;
           for(j=0;j<d[0];j++)
            prod *= ( x[i+n[0]*j] <= u[j]);
           M[i] = (double)prod;
           sum+= prod;
       }

     //  printf("\n sum = %d\n", sum);
       p = ((double)sum)/((double)n[0]);


     temp = 0.0;
     for(i=0;i<n[0];i++)
        {
          M[i] = M[i]-p; // centered now!
          temp += M[i];
          cumsum[i] = temp;
      }

   }

   void StatST(double *M, int *n, double *S, double *T)
{
    int i,j;
    double s1, s2, s3;
    double sn = 1.0/sqrt((double)n[0]);
    double n1 = 1.0/((double)n[0]);

    for(i=0;i<n[0];i++)
       {
           s1 = 0.0;
           s2 = 0.0;
           for(j=0;j<n[0];j++)
           {
             s3 = M[i+j*n[0]];
             s2 = maxi(s2,fabs(s3));
             s1 += s3*s3;
           }
           S[i] = s1*n1*n1;
           T[i] = s2*sn;
       }


}

 void cpCopulaStats(double *x, int *n, int *d, double *M, double *S, double *T)
   {
       int i,j,k;


       double *cum0 = calloc( n[0], sizeof(double));
       double *M0  = calloc( n[0], sizeof(double));
       double *u   = calloc(d[0], sizeof(double));
       double *cumsumcentered = calloc( n[0]*n[0], sizeof(double));
      // printf("n = %d d = %d\n", n[0],d[0]);



       for(j=0;j<n[0];j++)
       {
        for(k=0;k<d[0];k++)
            u[k] = x[j+n[0]*k];

         empcdf(x,n,d, u, M0, cum0);

           for(i=0;i<n[0];i++)
             {
               cumsumcentered[i+j*n[0]] = cum0[i] ;
               M[i+j*n[0]] = M0[i] ; // M is centered
             }

       }

        StatST(cumsumcentered,n,S,T);

       free(M0); free(cum0); free(cumsumcentered); free(u);
}


void cpCopulaStatsMult(double *M, double *xi, double *s, int *n, double *S, double *T)
   {
       int i,j;

       double s1;

      double *B = calloc( n[0]*n[0], sizeof(double));
      double *beta = calloc( n[0],  sizeof(double));


      for(j=0;j<n[0];j++)
       {
           s1 = 0.0;
          for(i=0;i<n[0];i++)
           {
             s1 += xi[i]*M[i+j*n[0]];
             beta[i]=s1;
           }

           for(i=0;i<n[0];i++)
             B[i+j*n[0]] = beta[i]-s[i]*beta[n[0]-1];


       }

       StatST(B,n,S,T);


       free(beta); free(B);
}


 void cpChangePointDStat(double *X, int *n, int *d, double *S, double *T)
    {
     /* Computes S and T statistics */

      double *U1, *U2;
      double *MM1, *MM2, *z1, *z2, *r1, *r2, *D;

      double s, temp,  B1, B2;

      int i1, i2, j, k, t;


   D = calloc(n[0], sizeof(double));
   S[0] = 0.0;
   T[0] = 0.0;


   for(i1=1;i1<n[0];i1++)
   {

    i2 = n[0]-i1;

    U1  = calloc(d[0]*i1 , sizeof(double));
    z1  = calloc(i1      , sizeof(double));
    r1  = calloc(i1      , sizeof(double));
    MM1 = calloc(i1      , sizeof(double));
    U2  = calloc(d[0]*i2 , sizeof(double));
    z2  = calloc(i2      , sizeof(double));
    r2  = calloc(i2      , sizeof(double));
    MM2 = calloc(i2      , sizeof(double));


      for(k=0;k<d[0];k++)
      {
         for(j=0;j<i1;j++)
             z1[j] = X[j+k*n[0]];

         rank(z1, r1,i1);

         for(j=0;j<i1;j++)
             U1[j+k*i1] = r1[j];

         for(j=0;j<i2;j++)
             z2[j] = X[i1+j+k*n[0]];

         rank(z2, r2,i2)  ;

         for(j=0;j<i2;j++)
             U2[j+k*i2] = r2[j];
      }



     for(t=0;t<n[0];t++)
     {
        for(j=0;j<i1;j++)
        {
          temp = 1.0;
          for(k=0;k<d[0];k++)
           temp *= (U1[j+k*i1] <= i1*X[t+k*n[0]]);

          MM1[j] = temp;
        }

        B1 = sum(MM1,i1)/sqrt((double)n[0]);



        /* tilde */
        for(j=0;j<i2;j++)
        {
          temp = 1.0;
          for(k=0;k<d[0];k++)
           temp *= (U2[j+k*i2] <= i2*X[t+k*n[0]]);

          MM2[j] = temp;
        }
         B2 = sum(MM2, i2)/sqrt((double)n[0]);




           D[t] = (i2*B1 - i1*B2)/( (double)n[0] );
     }
      T[i1] = maxvec(D,n[0]);
      s = 0.0;
      for(t=0;t<n[0];t++)
          s += D[t]*D[t];

      S[i1] = s/((double)n[0]);


   free(U1); free(U2); free(MM1); free(MM2); free(z1); free(z2); free(r1); free(r2);

   }

  free(D);

    }





void gradientCop(double *x, int *n, int *d, double *u, double *grad, double *MC, double *MC1)
   {
       int i,j,k;
       int pp, prod1, prod2, sum1, sum2 ,prod, sum, *s;
       double p, h, *v1,*v2;

       h = 1.0/sqrt((double)n[0]);

       v1 = calloc( d[0], sizeof(double));
       v2 = calloc( d[0], sizeof(double));
       s  = calloc( d[0], sizeof(int));

       for(j=0;j<d[0];j++)
         s[j] = 0;
      //printf("n = %d d = %d h = %f\n", n[0],d[0],h);


      for(k=0;k<d[0];k++)
       {

          sum1 = 0;
          sum2 = 0;
           for(j=0;j<d[0];j++)
           {
               v1[j]=u[j];
               v2[j]=u[j];
           }

           v1[k] = u[k] - h;
           v2[k] = u[k] + h;

         for(i=0;i<n[0];i++)
          {
           prod1 = 1;
           prod2 = 1;
           for(j=0;j<d[0];j++)
            {
                 prod1 *= ( x[i+n[0]*j] <= v1[j]);
                 prod2 *= ( x[i+n[0]*j] <= v2[j]);

            }
           sum1+= prod1;
           sum2+= prod2;
          }

       grad[k] = 0.5*(sum2-sum1)/((double)n[0])/h;

       }

       sum = 0;

       for(i=0;i<n[0];i++)
          {
           prod = 1;
           for(j=0;j<d[0];j++)
            {
               pp =  ( x[i+n[0]*j] <= u[j]);
               prod *= pp;
               MC1[i+n[0]*j] = pp;
               s[j] = s[j] + pp;
            }
           sum+= prod;

           MC[i] = (double)prod;
          }

        p = ((double)sum)/((double)n[0]);
        for(j=0;j<d[0];j++)
            v1[j] = ((double)s[j])/((double)n[0]);

      for(i=0;i<n[0];i++)
        {
            MC[i] = MC[i]-p;
            for(j=0;j<d[0];j++)
                MC1[i+n[0]*j] = MC1[i+n[0]*j]-v1[j];
        }


        free(v1); free(v2); free(s);
   }


   void cpCopulaStatsBucher(double *x, int *n, int *d, double *MC, double *MC1, double *grad)
   {
       int i,j,k;


       double *MC0, *MC10 , *u, *grad0;
      // printf("n = %d d = %d\n", n[0],d[0]);

      MC0    = calloc( n[0]       ,sizeof(double));
      MC10   = calloc( n[0]* d[0] ,sizeof(double));
      grad0  = calloc( d[0]       ,sizeof(double));
       u     = calloc(d[0]        ,sizeof(double));

       for(j=0;j<n[0];j++)
       {
          for(k=0;k<d[0];k++)
            u[k] = x[j+n[0]*k];

           gradientCop(x, n, d, u, grad0, MC0, MC10);

          for(k=0;k<d[0];k++)
            grad[k+j*d[0]] = grad0[k] ;

          for(i=0;i<n[0];i++)
            {
              MC[i+j*n[0]] = MC0[i] ;
              for(k=0;k<d[0];k++)
                MC1[i+j*n[0]*d[0]+k*n[0]] = MC10[i+n[0]*k] ;
            }

       }




       free(MC0); free(MC10); free(u); free(grad0);
}


void cpCopulaStatsMultBucherNonSeq(double *MC, double *MC1, double *grad, double *xi, double *s, int *n, int *d, double *S, double *T)
   {
       int i,j,k;

       double s1,s2;
       double *beta, *beta1, *B, *ss;


      B     = calloc( n[0]*n[0] ,sizeof(double));
      beta  = calloc( n[0]      ,sizeof(double));
      beta1 = calloc( n[0]*d[0] ,sizeof(double));
      ss    = calloc( d[0]      ,sizeof(double));


      for(j=0;j<n[0];j++)
       {
           s1 = 0.0;
          for(k=0;k<d[0];k++)
              ss[k] = 0.0;

          for(i=0;i<n[0];i++)
           {
             s1 += xi[i]*MC[i+j*n[0]];
             beta[i]=s1;

               for(k=0;k<d[0];k++)
               {
                 ss[k] = ss[k] + xi[i]*MC1[i + k*n[0] + j*n[0]*d[0]];
                 beta1[i+k*n[0]] = ss[k];

               }

           }


           for(i=0;i<n[0];i++)
             {
                 B[i+j*n[0]] = beta[i]-s[i]*beta[n[0]-1];

                 s2 = 0.0;
                 for(k=0;k<d[0];k++)
                    s2 += (beta1[i+k*n[0]]-s[i]*beta1[n[0]-1+k*n[0]] )*grad[k+j*d[0]];

                  B[i+j*n[0]] = B[i+j*n[0]] - s2;
               //

             }
       }


       StatST(B,n,S,T);


     free(beta); free(beta1); free(B); free(ss);

           //      printf("toto\n");
}

void cpCopulaStatsMultBucherSeq(double *U, double *grad, double *xi, int *n, int *d, double *S, double *T)
   {
       /* Implementation of the non-sequential method with multipliers */

      double *U1, *U2, *MM1, *MM2, *MM1m, *MM2m, *Dcheck, *z1, *z2, *r1, *r2;

      double s1,s2, temp, C1, C1m, C2, C2m, B1, B2, h;

      int i1, i2, j, k, t;


   Dcheck = calloc(n[0], sizeof(double));
   h = 1.0/sqrt((double)n[0]);

   for(i1=1;i1<n[0];i1++)
   {


     i2 = n[0]-i1;

    U1  = calloc(d[0]*i1 ,sizeof(double));
    z1  = calloc(i1      ,sizeof(double));
    r1  = calloc(i1      ,sizeof(double));
    MM1 = calloc(i1      ,sizeof(double));
    MM1m= calloc(i1      ,sizeof(double));
    U2  = calloc(d[0]*i2 ,sizeof(double));
    z2  = calloc(i2      ,sizeof(double));
    r2  = calloc(i2      ,sizeof(double));
    MM2 = calloc(i2      ,sizeof(double));
    MM2m= calloc(i2      ,sizeof(double));


      for(k=0;k<d[0];k++)
      {
         for(j=0;j<i1;j++)
             z1[j] = U[j+k*n[0]];

         rank(z1, r1,i1);

         for(j=0;j<i1;j++)
             U1[j+k*i1] = r1[j];

         for(j=0;j<i2;j++)
             z2[j] = U[i1+j+k*n[0]];

         rank(z2, r2,i2)  ;

         for(j=0;j<i2;j++)
             U2[j+k*i2] = r2[j];
      }

    for(t=0;t<n[0];t++)
     {
        for(j=0;j<i1;j++)
        {
          temp = 1.0;
          for(k=0;k<d[0];k++)
           temp *= (U1[j+k*i1] <= i1*U[t+k*n[0]]);

          MM1[j] = temp;
        }

        C1 = mean(MM1,i1);

        for(j=0;j<i1;j++)
          MM1[j] = xi[j]*(MM1[j]-C1);

        B1 = h*sum(MM1,i1);

       /******************************* 1-d processes */
        s1 = 0.0;
       for(k=0;k<d[0];k++)
         {
           for(j=0;j<i1;j++)
             MM1m[j] = (U1[j+k*i1] <= i1*U[t+k*n[0]]);


           C1m = mean(MM1m, i1);


           for(j=0;j<i1;j++)
             MM1m[j] = xi[j]*(MM1m[j]-C1m);

         s1 +=  sum(MM1m,i1)* grad[k+t*d[0]];
       }

        s1 = s1*h;

        /* tilde */
        for(j=0;j<i2;j++)
        {
          temp = 1.0;
          for(k=0;k<d[0];k++)
           temp *= (U2[j+k*i2] <= i2*U[t+k*n[0]]);

          MM2[j] = temp;
        }

        C2 = mean(MM2,i2);
         for(j=0;j<i2;j++)
           MM2[j] = xi[i1+j]*(MM2[j]-C2);

        B2 = h*sum(MM2, i2);

        /**************************** 1-d processes */
           s2 = 0.0;
        for(k=0;k<d[0];k++)
         {
           for(j=0;j<i2;j++)
             MM2m[j] = (U2[j+k*i2] <= i2*U[t+k*n[0]]);

           C2m = mean(MM2m, i2);


           for(j=0;j<i2;j++)
             MM2m[j] = xi[i1+j]*(MM2m[j]-C2m);

          s2 +=  sum(MM2m,i2)* grad[k+t*d[0]];
          }

           s2 = s2*h;
           Dcheck[t] = (i2*(B1-s1) - i1*(B2-s2))/( (double)n[0] );


     }



      T[i1] = maxvec(Dcheck,n[0]);
      s1 = 0.0;
      for(t=0;t<n[0];t++)
          s1 += Dcheck[t]*Dcheck[t];

      S[i1] = s1/((double)n[0]);

   free(U1); free(U2);
   free(MM1); free(MM2); free(MM1m); free(MM2m);
   free(z1); free(z2);free(r1); free(r2);


   }
  free(Dcheck);
   }
