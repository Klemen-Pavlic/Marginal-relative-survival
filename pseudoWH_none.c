/* 
	Author: Klemen Pavlic 
*/


void pseudoWH_none(double *Y, double *D, double *time, double *obsT, double *status, double *S_E, double *Sp, int *N, int *NT, int *neval, int *ind){
	int t, k, i_t, j;

	double S = 1; // KM on the whole sample
	
	double Sk[*N]; // KM on reduced sample
	
	double S_Oi[*N]; // pseudo observations

	double na_k;

	for (k=0; k<*N; k++){
		Sk[k] = 1;	
	}
	
	for (t=0;t<*neval;t++)
		{	
			i_t = ind[t];
			if (i_t == -1){
					//before the first event
					for (k=0;k<*N;k++){
							S_Oi[k] = 1;
							S_E[t] = S_E[t] + S_Oi[k] / Sp[t*(*N) + k];
						}
				}
			else { //after or at the first event
					if (t == 0){ // at the first event
							j = 0;						
							while (j < i_t){
									S = S * (1 - D[j] / Y[j]);
									for (k=0;k<*N;k++){
										/* compute KM estimators on the reduced samples */
										if (obsT[k]>time[j])
											{
												/* decrease the number at risk because individual k was at risk at time[t] */
												na_k = D[j] / (Y[j] - 1);
											}
										else{
												if (obsT[k]==time[j])
													{
														/* decrease the number of events if k was an event, and decrease the number at risk because
														k was in the risk set at time[t] */
														na_k = (D[j] - status[k]) / (Y[j] - 1);  
													}
												else{ // obsT[k]<time[j]
														/* do nothing since k is no longer in the risk set */
														na_k = D[j] / Y[j];
													}
											}
										/* Compute the estimates */
										Sk[k] = Sk[k] * (1 - na_k);
										S_Oi[k] = (*N) * S - (*N-1) * Sk[k];	
										}
									j++;
								}
							/*j==i_t*/
							S = S * (1 - D[i_t] / Y[i_t]);
							for (k=0;k<*N;k++){
								/* compute KM estimators on the reduced samples */
								if (obsT[k]>time[i_t])
									{
										/* decrease the number at risk because individual k was at risk at time[t] */
										na_k = D[i_t] / (Y[i_t] - 1);
									}
								else{
										if (obsT[k]==time[i_t])
											{
												/* decrease the number of events if k was an event, and decrease the number at risk because
												k was in the risk set at time[t] */
												na_k = (D[i_t] - status[k]) / (Y[i_t] - 1);
											}
										else{ // obsT[k]<time[i_t]
												/* do nothing since k is no longer in the risk set */
												na_k = D[i_t] / Y[i_t];
											}
									}
								/* Compute the estimates */
								Sk[k] = Sk[k] * (1 - na_k);
								S_Oi[k] = (*N) * S - (*N-1) * Sk[k];
								S_E[t] = S_E[t] + S_Oi[k]/Sp[t*(*N) + k];
								}
						}
					else {
							if (i_t == ind[t-1]){ /*between two event times*/
									for (k=0;k<*N;k++){
										/* Compute the estimates */
										/*S_Oi[k] = (*N) * S - (*N-1) * Sk[k];	*/
										S_E[t] = S_E[t] + S_Oi[k]/Sp[t*(*N) + k];
										}
								}
							else if (i_t - ind[t-1] > 1) {
									j = ind[t-1] + 1;						
									while (j < i_t){
											S = S * (1 - D[j] / Y[j]);
											for (k=0;k<*N;k++){
												/* compute KM estimators on the reduced samples */
												if (obsT[k]>time[j])
													{
														/* decrease the number at risk because individual k was at risk at time[t] */
														na_k = D[j] / (Y[j] - 1);
													}
												else{
														if (obsT[k]==time[j])
															{
																/* decrease the number of events if k was an event, and decrease the number at risk because
																k was in the risk set at time[t] */
																na_k = (D[j] - status[k]) / (Y[j] - 1);  
															}
														else{ // obsT[k]<time[j]
																/* do nothing since k is no longer in the risk set */
																na_k = D[j] / Y[j];
															}
													}
												/* Compute the estimates */
												Sk[k] = Sk[k] * (1 - na_k);
												S_Oi[k] = (*N) * S - (*N-1) * Sk[k];	
												}
											j++;
										}
									/*j==i_t*/
									S = S * (1 - D[j] / Y[j]);
									for (k=0;k<*N;k++){
										/* compute KM estimators on the reduced samples */
										if (obsT[k]>time[i_t])
											{
												/* decrease the number at risk because individual k was at risk at time[t] */
												na_k = D[i_t] / (Y[i_t] - 1);
											}
										else{
												if (obsT[k]==time[i_t])
													{
														/* decrease the number of events if k was an event, and decrease the number at risk because
														k was in the risk set at time[t] */
														na_k = (D[i_t] - status[k]) / (Y[i_t] - 1);
													}
												else{ // obsT[k]<time[i_t]
														/* do nothing since k is no longer in the risk set */
														na_k = D[i_t] / Y[i_t];
													}
											}
										/* Compute the estimates */
										Sk[k] = Sk[k] * (1 - na_k);
										S_Oi[k] = (*N) * S - (*N-1) * Sk[k];	
										S_E[t] = S_E[t] + S_Oi[k]/Sp[t*(*N) + k];
										}
								}						
							else { /* at event times*/
									S = S * (1 - D[i_t] / Y[i_t]);
									for (k=0;k<*N;k++){
										/* compute KM estimators on the reduced samples */
										if (obsT[k]>time[i_t])
											{
												/* decrease the number at risk because individual k was at risk at time[t] */
												na_k = D[i_t] / (Y[i_t] - 1);
											}
										else{
												if (obsT[k]==time[i_t])
													{
														/* decrease the number of events if k was an event, and decrease the number at risk because
														k was in the risk set at time[t] */
														na_k = (D[i_t] - status[k]) / (Y[i_t] - 1);
													}
												else{ // obsT[k]<time[i_t]
														/* do nothing since k is no longer in the risk set */
														na_k = D[i_t] / Y[i_t];
													}
											}
										/* Compute the estimates */
										Sk[k] = Sk[k] * (1 - na_k);
										S_Oi[k] = (*N) * S - (*N-1) * Sk[k];	
										S_E[t] = S_E[t] + S_Oi[k]/Sp[t*(*N) + k];
										}
								}
						}
				}
			S_E[t] = S_E[t] / (*N);
		}
	}
