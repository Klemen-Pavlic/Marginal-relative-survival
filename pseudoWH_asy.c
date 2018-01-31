/* 
	Author: Klemen Pavlic 
*/


void pseudoWH_asy(double *Y, double *D, double *time, double *obsT, double *status, double *S, double *varS, 
														double *Sk, double *varSk, double *covSSk, int *N, int *NT)
{
	int t, k;
	S[0] = 1;
	varS[0] = 0;
	double na_k, v_k_int, c_k_int;
	for (t=0;t<*NT;t++)
		{	
			// KM on the whole sample
			S[t+1] = S[t] * (1 - D[t]/Y[t]);
			// var of KM on the whole sample
			varS[t+1] = varS[t] + D[t] / ( (Y[t] - D[t]) * Y[t] );
			for (k=0;k<*N;k++)
				{	
					if (t==0)
						{
							/* set the values at time zero */
							Sk[k] = 1;
							varSk[k] = 0;
							covSSk[k] = 0;
						}
					/* compute KM estimators on the reduced samples, their variances and covariances with the KM on the whole sample */
					if (obsT[k]>time[t])
						{
							/* decrease the number at risk because individual k was at risk at time[t] */
							na_k = D[t] / (Y[t] - 1);
							v_k_int = D[t] / ( (Y[t] - 1 - D[t]) * (Y[t] - 1) );
							c_k_int = (Y[t] - 1) * D[t] / ( Y[t] * Y[t] * ( Y[t] - 1 - D[t] ) );
						}
					else{
							if (obsT[k]==time[t])
								{
									/* decrease the number of events if k was an event, and decrease the number at risk because
									k was in the risk set at time[t] */
									na_k = (D[t] - status[k]) / (Y[t] - 1);
									v_k_int = (D[t] - status[k]) / ( ( Y[t] - 1 - (D[t]-status[k]) ) * (Y[t] - 1) );
			            c_k_int = (Y[t] - 1) * D[t] / ( Y[t] * Y[t] * ( Y[t] - 1 - ( D[t] - status[k] ) ) );  
								}
							else{ // obsT[k]<time[t]
									/* do nothing since k is no longer in the risk set */
									na_k = D[t] / Y[t];
									v_k_int = D[t] / ( (Y[t] - D[t]) * Y[t] );
            			c_k_int = D[t] / ( Y[t] * ( Y[t] - D[t] )); // Y[t] * D[t] / ( Y[t] * Y[t] * ( Y[t] - D[t] ) ) reduces to
								}
						}
					/* Compute the estimates */
					Sk[k+(*N)*(t+1)] = Sk[k+(*N)*t] * (1 - na_k);
					varSk[k+(*N)*(t+1)] = varSk[k+(*N)*t] + v_k_int;
					covSSk[k+(*N)*(t+1)] = covSSk[k+(*N)*t] + c_k_int;
				}
		}
}
/*
void pseudoWH_pre(double *Y, double *D, double *time, double *obsT, double *status, double *S, double *varS,
														double *Sk, double *covSkSl, double *covSSk,int *N, int *NT)
{

}
*/

#include <stdio.h>
int main()
{	
	double time[] = { 71, 972, 1032, 1269, 1394, 1398, 1654};
	double Y[] = { 10, 9, 8, 7, 6, 5, 4};
	double D[] = { 1, 1, 1, 1, 1, 1, 1};
	int NT = 7;
	double obsT[] = { 71, 972, 1032, 1269, 1394, 1398, 1654, 6275, 6574, 7682};
	double status[] = { 1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
	int N = 10;
	double S[NT+1];
	double varS[NT+1];	
	double Sk[N*(NT+1)];
	double varSk[N*(NT+1)];
	double covSSk[N*(NT+1)];

	double *time_p = &time[0];
	double *Y_p = &Y[0];
	double *D_p = &D[0]; 
	int *NT_p;
	NT_p = &NT;
	double *obsT_p = &obsT[0];
	double *status_p = &status[0]; 
	int *N_p;
	N_p = &N;
	double *S_p = &S[0];
	double *varS_p = &varS[0];
	double *Sk_p = &Sk[0];
	double *varSk_p = &varSk[0];
	double *covSSk_p = &covSSk[0];	
	
	int l = sizeof(S)/sizeof(S_p);
	int l2 = sizeof(Sk)/sizeof(Sk_p);
	int i;

	pseudoWH_asy(Y_p, D_p, time_p, obsT_p, status_p, S_p, varS_p, Sk_p, varSk_p, covSSk_p, N_p, NT_p);
	
	
	printf("\nPo evaluaciji.\n");	
	printf("Vrednost S je:\n");
	for (i=0;i<l;i++){
		printf("%f ",S[i]);
	}
	printf("\nVrednost varS je:\n");
	for(i=0;i<l;i++){
		printf("%f ",varS[i]);		
		//printf("%f\n", *(varS_p+i));
	}
	printf("\nVrednost Sk je:\n");
	for (i=0;i<l2;i++){
		printf("%f ", Sk[i]);
	}
	printf("\nVrednost varSk je:\n");
	for(i=0;i<l2;i++){
		printf("%f ",*(varSk_p + i));
	}
	printf("\nVrednost covSSk je:\n");
	for (i=0;i<l2;i++){
		printf("%f ", *(covSSk_p + i));
		//printf("%f\n", covSSk[i]);
	}
	
	return 0;
}

