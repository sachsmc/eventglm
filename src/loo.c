#include <math.h>
#include <R.h>

void loo_surv2(double *Y,
	      double *D,
	      double *time,
	      double *obsT,
	      double *status,
	      double *S,
	      int *N,
	      int *NT,
	      int *Tdex){
  int k, t;
  double na,pl;
  for (k=0; k<*N;k++){
    /* Rprintf("\n"); */
    /* compute the Nelson-Aalen estimate */
    pl=1;
    for (t=0; t<*NT;t++){
      if (obsT[k]>time[t]){
	/* decrease the number at risk
	   because individual k was at risk
	   at time[t]
	 */
	if(Y[t] == 1){
	  na = D[t];
	} else {
	  na = D[t]/(Y[t]-1);
	}
      }
      else{
	if (obsT[k]==time[t]){
	  /*
	    decrease the number of events
	    if k was an event,
	    and decrease the number at risk
	    because k was in the risk set at
	    time[t]
	  */
	  if(Y[t] == 1) {
	    na = (D[t]-status[k]);
	  } else {
	    na = (D[t]-status[k])/(Y[t]-1);
	  }

	}
	else{
	  /* do nothing */
	  na = D[t]/Y[t];
	}
      }
      /* compute the product-limit estimate */
      pl *= (1-na);
      if(t == *Tdex) {
      S[k]=pl;
      }
      /* Rprintf("t=%d\tk=%d\tD[t]=%1.2f\tY[t]=%1.2f\tna=%1.2f\tS[k](t)=%1.2f\n",t,k,D[t],Y[t],na,S[k+(*N)*t]); */
    }
  }
}

void loo_surv(double *Y,
               double *D,
               double *time,
               double *obsT,
               double *status,
               double *S,
               int *N,
               int *NT,
               int *Tdex){
    int k, t;
    double na,pl;
    for (k=0; k<*N;k++){
        /* Rprintf("\n"); */
        /* compute the Nelson-Aalen estimate */
        pl=1;
        for (t=0; t<=*Tdex;t++){
            if (obsT[k]>time[t]){
                /* decrease the number at risk
                 because individual k was at risk
                 at time[t]
                 */
                if(Y[t] == 1){
                  na = D[t];
                } else {
                  na = D[t]/(Y[t]-1);
                }
            }
            else{
                if (obsT[k]==time[t]){
                    /*
                     decrease the number of events
                     if k was an event,
                     and decrease the number at risk
                     because k was in the risk set at
                     time[t]
                     */
                    if(Y[t] == 1){
                      na = (D[t]-status[k]);
                    } else {
                      na = (D[t]-status[k])/(Y[t]-1);
                    }
                }
                else{
                    /* do nothing */
                    na = D[t]/Y[t];
                }
            }
            /* compute the product-limit estimate */
            pl *= (1-na);

            S[k+(*N)*t]=pl;

            /* Rprintf("t=%d\tk=%d\tD[t]=%1.2f\tY[t]=%1.2f\tna=%1.2f\tS[k](t)=%1.2f\n",t,k,D[t],Y[t],na,S[k+(*N)*t]); */
        }
    }
}


void loo_comprisk2(double *Y,
                  double *D,
                  double *D0,
                  double *time,
                  double *obsT,
                  double *status,
                  double *status0,
                  double *F,
                  int *N,
                  int *NT,
                  int *Tdex){
    int k, t;
    double na,aj,na0,ls;
    for (k=0; k<*N;k++){
        /* compute the Nelson-Aalen estimate */
        aj=0;
        ls=1;
        for (t=0; t<*NT;t++){
            if(t > 0) {
                if(obsT[k] > time[t - 1]) {
                  if(Y[t - 1] == 1) {
                    na0 = D0[t - 1];
                  } else {
                    na0 = D0[t - 1]/(Y[t - 1]-1);
                  }
                }
                else {
                    if(obsT[k]==time[t-1]) {
                      if(Y[t - 1]== 1) {
                        na0 = (D0[t-1]-status0[k]);
                      } else {
                        na0 = (D0[t-1]-status0[k])/(Y[t-1]-1);
                      }
                    }
                    else {
                        na0 = D0[t-1]/Y[t-1];
                    }

                }

            } else {
                na0 = 0;
            }


            if (obsT[k]>time[t]){
                /* decrease the number at risk
                 because k was in the risk set at time[t]
                 */

              if(Y[t] == 1) {
                na = D[t];
              } else {
                na = D[t]/(Y[t]-1);
              }



            }
            else{
                if (obsT[k]==time[t]){
                    /*
                     decrease the number of events
                     if k was an event,
                     and decrease the number at risk
                     because k was in the risk set at
                     time[t]
                     */

                  if(Y[t] == 1) {
                    na = D[t] - status[k];
                  } else {
                    na = (D[t]-status[k])/(Y[t]-1);
                  }

                }
                else{
                    /* do nothing */
                    na = D[t]/Y[t];

                }
            }
            ls *= (1-na0);
            /* compute the Aalen-Johansen estimate */
            aj += ls*na;

            if(t == *Tdex){
            F[k]=aj;
            }
        }
    }
}


void loo_comprisk(double *Y,
                   double *D,
                   double *D0,
                   double *time,
                   double *obsT,
                   double *status,
                   double *status0,
                   double *F,
                   int *N,
                   int *NT,
                   int *Tdex){
    int k, t;
    double na,aj,na0,ls;
    for (k=0; k<*N;k++){
        /* compute the Nelson-Aalen estimate */
        aj=0;
        ls=1;
        for (t=0; t<=*Tdex;t++){
            if(t > 0) {
                if(obsT[k] > time[t - 1]) {
                  if(Y[t - 1] == 1) {
                    na0 = D0[t - 1];
                  } else {
                    na0 = D0[t - 1]/(Y[t - 1]-1);
                  }
                }
                else {
                    if(obsT[k]==time[t-1]) {

                      if(Y[t-1] == 1) {
                        na0 = (D0[t-1]-status0[k]);
                      } else {
                        na0 = (D0[t-1]-status0[k])/(Y[t-1]-1);
                      }

                    }
                    else {
                        na0 = D0[t-1]/Y[t-1];
                    }

                }

            } else {
                na0 = 0;
            }


            if (obsT[k]>time[t]){
                /* decrease the number at risk
                 because k was in the risk set at time[t]
                 */
                if(Y[t] == 1) {
                  na = D[t];
                } else {
                  na = D[t]/(Y[t]-1);
                }




            }
            else{
                if (obsT[k]==time[t]){
                    /*
                     decrease the number of events
                     if k was an event,
                     and decrease the number at risk
                     because k was in the risk set at
                     time[t]
                     */
                    if(Y[t] == 1) {
                      na = D[t] - status[k];
                    } else {
                      na = (D[t]-status[k])/(Y[t]-1);
                    }


                }
                else{
                    /* do nothing */
                    na = D[t]/Y[t];

                }
            }
            ls *= (1-na0);
            /* compute the Aalen-Johansen estimate */
            aj += ls*na;

            F[k+(*N)*t]=aj;

        }
    }
}
