/* 
	Author: Klemen Pavlic 
*/

void pseudoWH_pre2(double *Y, double *D, double *time, double *obsT, double *status, double *S_E, double *varS_E, double *Sp, int *N, int *NT, int *neval, int *ind, int *varSkMethod, double *covSkSl){
	int t, k, l, i_t, j;

	double S = 1; // KM on the whole sample
	double varS = 0; // the integral part of the variance
	
	double Sk[*N]; // KM on reduced sample

	double covSSk[*N]; // the integral part of the covariance
	//double covSkSl[(*N) * (*N+1) / 2]; // the integral part of the covariance 
	double S_Oi[*N]; // pseudo observations

	double v1, v2, v3; 
	double na_k, c_k_int, c_kl_int;
	
	for (k=0;k<*N;k++){
		Sk[k] = 1;
		//covSkSl[k] = 0;
		covSSk[k] = 0;		
	}
	//for (k=*N;k<(*N) * (*N+1) / 2;k++){covSkSl[k] = 0;}	

	for (t=0;t<*neval;t++)
		{	
			v1 = 0;
			v2 = 0;
			v3 = 0;
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
									varS = varS + D[j] / ( (Y[j] - D[j]) * Y[j] );
									for (k=0;k<*N;k++){
										//v1 = v1 + 1 / Sp[t*(*N) + k];
										if (k==0) {
											if (obsT[k]>time[j])
												{
													// decrease the number at risk because individual k was at risk at time[j] 
													na_k = D[j] / (Y[j] - 1);
													c_k_int = ( (Y[j] - 1) * D[j] ) / ( Y[j] * Y[j] * ( Y[j] - 1 - D[j] ) );
												}
											else{
													if (obsT[k]==time[j])
														{
															// decrease the number of events if k was an event, and decrease the number at risk because k was in the risk set at time[j]
															na_k = (D[j] - status[k]) / (Y[j] - 1);
															c_k_int = ( (Y[j] - 1) * D[j] ) / ( Y[j] * Y[j] * ( Y[j] - 1 - ( D[j] - status[k] ) ) );  
														}
													else{ // obsT[k]<time[i_t]
															// do nothing since k is no longer in the risk set 
															na_k = D[j] / Y[j];
															c_k_int = D[j] / ( Y[j] * ( Y[j] - D[j] ) ); // ( Y[j] * D[j]) / ( Y[j] * Y[j] * ( Y[j] - D[j] ) ) reduces to
														}
												}
											Sk[k] = Sk[k] * (1 - na_k);
											covSSk[k] = covSSk[k] + c_k_int;
											} // end of k==0 case
										// l == k case
										if (*varSkMethod == 1)
											{//precise
												if (obsT[k] > time[j]){
					 	              c_kl_int = (Y[j] - D[j]) * (Y[j] - 1) * D[j] / ( ( (Y[j] - 1) - D[j] ) * ( (Y[j] - 1) - D[j] ) * Y[j] * Y[j] );
		      					    }
		        				    else {
		        				      if (obsT[k] == time[j]){
		        					      c_kl_int = (Y[j] - D[j]) * (Y[j] - 1) * D[j] / ( ( (Y[j] - 1) - (D[j] - status[k]) ) * ( (Y[j] - 1) - (D[j] - status[k]) ) * Y[j] * Y[j] );
		        					    }
		        				      else { // obsT[k] < time[j]
		        					      c_kl_int = D[j] / ( (Y[j] - D[j]) * Y[j]); // (Y[j] - D[j]) * Y[j] * D[j] / ( (Y[j] - D[j]) * (Y[j] - D[j]) * Y[j] * Y[j] ) reduces to 
		        					    }
		            				}
											}
										else 		
											{//Reduced sample KM
												if (obsT[k] > time[j]){
		              				c_kl_int = D[j] / ( (Y[j] - 1 - D[j]) * (Y[j] - 1) );
		            				}
		            				else {
		              				if (obsT[k] == time[j]){
		                				c_kl_int = (D[j] - status[k]) / ( ( Y[j] - 1 - (D[j] - status[k]) ) * (Y[j] - 1) );
		              				}
		              				else { // obsT[k] < time[j]
		                				c_kl_int = D[j] / ( (Y[j] - D[j]) * Y[j]);
		              				}
		            				}
											}
										covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] = covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] + c_kl_int;
										//v2 = v2 + Sk[k] * Sk[k] * covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] / Sp[t*(*N) + k] / Sp[t*(*N) + k];
										// l > k case
										for (l=k+1;l<*N;l++){ // beginning of l-loop
											if (k==0){
												if (obsT[l]>time[j])
													{
														// decrease the number at risk because individual l was at risk at time[j]
														na_k = D[j] / (Y[j] - 1);
														c_k_int = ( (Y[j] - 1) * D[j] ) / ( Y[j] * Y[j] * ( Y[j] - 1 - D[j] ) );
													}
												else{
														if (obsT[l]==time[j])
															{
																// decrease the number of events if l was an event, and decrease the number at risk because l was in the risk set at time[j] 
																na_k = (D[j] - status[l]) / (Y[j] - 1);
																c_k_int = ( (Y[j] - 1) * D[j] ) / ( Y[j] * Y[j] * ( Y[j] - 1 - ( D[j] - status[l] ) ) );  
															}
														else{ // obsT[l]<time[j]
																// do nothing since l is no longer in the risk set 
																na_k = D[j] / Y[j];
																c_k_int = D[j] / ( Y[j] * ( Y[j] - D[j] ) ); // ( Y[j] * D[j]) / ( Y[j] * Y[j] * ( Y[j] - D[j] ) ) reduces to
															}
													}
												Sk[l] = Sk[l] * (1 - na_k);
												covSSk[l] = covSSk[l] + c_k_int;
												}
											// k!= 0
											if (obsT[k] > time[j]){
								        if (obsT[l] > time[j]){
						          		c_kl_int = ( (Y[j] - D[j]) * (Y[j] - 2) * D[j] ) / ( (Y[j] - 1 - D[j]) * (Y[j] - 1 - D[j]) * Y[j] * Y[j] );
						        		}
						        		else {
						          		if (obsT[l] == time[j]){
						          		  c_kl_int = ( (Y[j] - D[j]) * (Y[j] - 2) * D[j] ) / ( (Y[j] - 1 - D[j]) * (Y[j] - 1 - (D[j] - status[l])) * Y[j] * Y[j] );
						          		}
						          		else { // obsT[l] < time[j]
						          		  c_kl_int = ( (Y[j] - 1) * D[j] ) / ( (Y[j] - 1 - D[j]) * Y[j] * Y[j] ); // ( (Y[j] - D[j]) * (Y[j] - 1) * D[j] ) / ( (Y[j] - 1 - D[j]) * (Y[j] - D[j]) * Y[j] * Y[j] ) reduces to 
						          		}
						        		}
						      		}
						      		else {
						        		if (obsT[k] == time[j]){
						          		if (obsT[l] > time[j]){
						            		c_kl_int = ( (Y[j] - D[j]) * (Y[j] - 2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k])) * (Y[j] - 1 - D[j]) * Y[j] * Y[j] );
						          		}
						          		else {
						            		if (obsT[l] == time[j]){
						              		c_kl_int = ( (Y[j] - D[j]) * (Y[j] - 2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k])) * (Y[j] - 1 - (D[j] - status[l])) * Y[j] * Y[j] );
						            		}
						            		else { // obsT[l] < time[j]
						              		c_kl_int = ( (Y[j] - 1) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k])) * Y[j] * Y[j] ); // ( (Y[j] - D[j]) * (Y[j] - 1) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k])) * (Y[j] - D[j]) * Y[j] * Y[j] ) reduces to 
						            		}
			 		                }	
								        }
								        else { // obsT[k] < time[j]
								          if (obsT[l] > time[j]){
								            c_kl_int = ( (Y[j] - 1) * D[j] ) / ( (Y[j] - 1 - D[j]) * Y[j] * Y[j] ); // ( (Y[j] - D[j]) * (Y[j] - 1) * D[j] ) / ( (Y[j] - D[j]) * (Y[j] - 1 - D[j]) * Y[j] * Y[j] ) reduces to 
								          }
								          else {
								            if (obsT[l] == time[j]){
								              c_kl_int = ( (Y[j] - 1) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[l])) * Y[j] * Y[j] ); // ( (Y[j] - D[j]) * (Y[j] - 1) * D[j] ) / ( (Y[j] - D[j]) * (Y[j] - 1 - (D[j] - status[l])) * Y[j] * Y[j] ) reduces to 
								            }
								            else { // obsT[l] < time[j]
								              c_kl_int = D[j] / ( (Y[j] - D[j]) * Y[j] ); // ( (Y[j] - D[j]) * Y[j] * D[j] ) / ( (Y[j] - D[j]) * (Y[j] - D[j]) * Y[j] * Y[j] ) reduces to 
								            }
								          }
								        }
								      }
											covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] = covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] + c_kl_int;
											//v2 = v2 + 2 * Sk[k] * Sk[l] * covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] / Sp[t*(*N) + k] / Sp[t*(*N) + l];
											} // end of l-loop
										// Compute the estimates 
										//S_Oi[k] = (*N) * S - (*N-1) * Sk[k];
										//S_E[t] = S_E[t] + S_Oi[k]/Sp[t*(*N) + k]; TUKAJ
										//v3 = v3 + S * Sk[k] * covSSk[k] / Sp[t*(*N) + k];
										} // end of k-loop
									j++;
								} // end of while loop over j
							//j==i_t
							S = S * (1 - D[i_t] / Y[i_t]);
							varS = varS + D[i_t] / ( (Y[i_t] - D[i_t]) * Y[i_t] );
							for (k=0;k<*N;k++){
								v1 = v1 + 1 / Sp[t*(*N) + k];
								if (k==0) {
									if (obsT[k]>time[i_t])
										{
											// decrease the number at risk because individual k was at risk at time[i_t]
											na_k = D[i_t] / (Y[i_t] - 1);
											c_k_int = ( (Y[i_t] - 1) * D[i_t] ) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - 1 - D[i_t] ) );
										}
									else{
											if (obsT[k]==time[i_t])
												{
													// decrease the number of events if k was an event, and decrease the number at risk because k was in the risk set at time[i_t] 
													na_k = (D[i_t] - status[k]) / (Y[i_t] - 1);
													c_k_int = ( (Y[i_t] - 1) * D[i_t] ) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - 1 - ( D[i_t] - status[k] ) ) );  
												}
											else{ // obsT[k]<time[i_t]
													// do nothing since k is no longer in the risk set 
													na_k = D[i_t] / Y[i_t];
													c_k_int = D[i_t] / ( Y[i_t] * ( Y[i_t] - D[i_t] ) ); // ( Y[i_t] * D[i_t]) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - D[i_t] ) ) reduces to
												}
										}
									Sk[k] = Sk[k] * (1 - na_k);
									covSSk[k] = covSSk[k] + c_k_int;
									} // end of k==0 case
								// l == k case
								if (*varSkMethod == 1)
									{//precise
										if (obsT[k] > time[i_t]){
			 	              c_kl_int = (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] / ( ( (Y[i_t] - 1) - D[i_t] ) * ( (Y[i_t] - 1) - D[i_t] ) * Y[i_t] * Y[i_t] );
      					    }
        				    else {
        				      if (obsT[k] == time[i_t]){
        					      c_kl_int = (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] / ( ( (Y[i_t] - 1) - (D[i_t] - status[k]) ) * ( (Y[i_t] - 1) - (D[i_t] - status[k]) ) * Y[i_t] * Y[i_t] );
        					    }
        				      else { // obsT[k] < time[i_t]
        					      c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t]) * Y[i_t]); // (Y[i_t] - D[i_t]) * Y[i_t] * D[i_t] / ( (Y[i_t] - D[i_t]) * (Y[i_t] - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to 
        					    }
            				}
									}
								else 		
									{//Reduced sample KM
										if (obsT[k] > time[i_t]){
              				c_kl_int = D[i_t] / ( (Y[i_t] - 1 - D[i_t]) * (Y[i_t] - 1) );
            				}
            				else {
              				if (obsT[k] == time[i_t]){
                				c_kl_int = (D[i_t] - status[k]) / ( ( Y[i_t] - 1 - (D[i_t] - status[k]) ) * (Y[i_t] - 1) );
              				}
              				else { // obsT[k] < time[i_t]
                				c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t]) * Y[i_t]);
              				}
            				}
									}
								covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] = covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] + c_kl_int;
								v2 = v2 + Sk[k] * Sk[k] * covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] / Sp[t*(*N) + k] / Sp[t*(*N) + k];
								// l > k case
								for (l=k+1;l<*N;l++){ // beginning of l-loop
									if (k==0){
										if (obsT[l]>time[i_t])
											{
												// decrease the number at risk because individual l was at risk at time[i_t] 
												na_k = D[i_t] / (Y[i_t] - 1);
												c_k_int = ( (Y[i_t] - 1) * D[i_t] ) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - 1 - D[i_t] ) );
											}
										else{
												if (obsT[l]==time[i_t])
													{
														// decrease the number of events if l was an event, and decrease the number at risk because l was in the risk set at time[i_t] 
														na_k = (D[i_t] - status[l]) / (Y[i_t] - 1);
														c_k_int = ( (Y[i_t] - 1) * D[i_t] ) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - 1 - ( D[i_t] - status[l] ) ) );  
													}
												else{ // obsT[l]<time[i_t]
														// do nothing since l is no longer in the risk set 
														na_k = D[i_t] / Y[i_t];
														c_k_int = D[i_t] / ( Y[i_t] * ( Y[i_t] - D[i_t] ) ); // ( Y[i_t] * D[i_t]) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - D[i_t] ) ) reduces to
													}
											}
										Sk[l] = Sk[l] * (1 - na_k);
										covSSk[l] = covSSk[l] + c_k_int;
										}
									// k!= 0
									if (obsT[k] > time[i_t]){
						        if (obsT[l] > time[i_t]){
				          		c_kl_int = ( (Y[i_t] - D[i_t]) * (Y[i_t] - 2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] );
				        		}
				        		else {
				          		if (obsT[l] == time[i_t]){
				          		  c_kl_int = ( (Y[i_t] - D[i_t]) * (Y[i_t] - 2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * (Y[i_t] - 1 - (D[i_t] - status[l])) * Y[i_t] * Y[i_t] );
				          		}
				          		else { // obsT[l] < time[i_t]
				          		  c_kl_int = ( (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * (Y[i_t] - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to 
				          		}
				        		}
				      		}
				      		else {
				        		if (obsT[k] == time[i_t]){
				          		if (obsT[l] > time[i_t]){
				            		c_kl_int = ( (Y[i_t] - D[i_t]) * (Y[i_t] - 2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k])) * (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] );
				          		}
				          		else {
				            		if (obsT[l] == time[i_t]){
				              		c_kl_int = ( (Y[i_t] - D[i_t]) * (Y[i_t] - 2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k])) * (Y[i_t] - 1 - (D[i_t] - status[l])) * Y[i_t] * Y[i_t] );
				            		}
				            		else { // obsT[l] < time[i_t]
				              		c_kl_int = ( (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k])) * Y[i_t] * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k])) * (Y[i_t] - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to 
				            		}
	 		                }	
						        }
						        else { // obsT[k] < time[i_t]
						          if (obsT[l] > time[i_t]){
						            c_kl_int = ( (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to 
						          }
						          else {
						            if (obsT[l] == time[i_t]){
						              c_kl_int = ( (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[l])) * Y[i_t] * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1 - (D[i_t] - status[l])) * Y[i_t] * Y[i_t] ) reduces to 
						            }
						            else { // obsT[l] < time[i_t]
						              c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t]) * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * Y[i_t] * D[i_t] ) / ( (Y[i_t] - D[i_t]) * (Y[i_t] - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to
						            }
						          }
						        }
						      }
									covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] = covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] + c_kl_int;
									v2 = v2 + 2 * Sk[k] * Sk[l] * covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] / Sp[t*(*N) + k] / Sp[t*(*N) + l];
									} // end of l-loop
								// Compute the estimates 
								S_Oi[k] = (*N) * S - (*N-1) * Sk[k];
								S_E[t] = S_E[t] + S_Oi[k]/Sp[t*(*N) + k];
								v3 = v3 + S * Sk[k] * covSSk[k] / Sp[t*(*N) + k];
								} // end of k-loop
						}
					else {
							if (i_t == ind[t-1]){ // between two event times
									for (k=0;k<*N;k++){
										v1 = v1 + 1 / Sp[t*(*N) + k];
										v3 = v3 + S * Sk[k] * covSSk[k] / Sp[t*(*N) + k];
										v2 = v2 + Sk[k] * Sk[k] * covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] / Sp[t*(*N) + k] / Sp[t*(*N) + k];
										for (l=k+1;l<*N;l++){
											v2 = v2 + 2 * Sk[k] * Sk[l] * covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] / Sp[t*(*N) + k] / Sp[t*(*N) + l];
											}
										// Compute the estimates 
										// S_Oi[k] = (*N) * S - (*N-1) * Sk[k];
										S_E[t] = S_E[t] + S_Oi[k]/Sp[t*(*N) + k];
										}
								}
							else if (i_t - ind[t-1] > 1) {
									j = ind[t-1] + 1;						
									while (j < i_t){
											S = S * (1 - D[j] / Y[j]);
											varS = varS + D[j] / ( (Y[j] - D[j]) * Y[j] );
											for (k=0;k<*N;k++){
											//v1 = v1 + 1 / Sp[t*(*N) + k];
											if (k==0) {
												if (obsT[k]>time[j])
													{
														// decrease the number at risk because individual k was at risk at time[j] 
														na_k = D[j] / (Y[j] - 1);
														c_k_int = ( (Y[j] - 1) * D[j] ) / ( Y[j] * Y[j] * ( Y[j] - 1 - D[j] ) );
													}
												else{
														if (obsT[k]==time[j])
															{
																// decrease the number of events if k was an event, and decrease the number at risk because k was in the risk set at time[j] 
																na_k = (D[j] - status[k]) / (Y[j] - 1);
																c_k_int = ( (Y[j] - 1) * D[j] ) / ( Y[j] * Y[j] * ( Y[j] - 1 - ( D[j] - status[k] ) ) );  
															}
														else{ // obsT[k]<time[j]
																// do nothing since k is no longer in the risk set 
																na_k = D[j] / Y[j];
																c_k_int = D[j] / ( Y[j] * ( Y[j] - D[j] ) ); // ( Y[j] * D[j]) / ( Y[j] * Y[j] * ( Y[j] - D[j] ) ) reduces to
															}
													}
												Sk[k] = Sk[k] * (1 - na_k);
												covSSk[k] = covSSk[k] + c_k_int;
												} // end of k==0 case
											// l == k case
											if (*varSkMethod == 1)
												{//precise
													if (obsT[k] > time[j]){
						 	              c_kl_int = (Y[j] - D[j]) * (Y[j] - 1) * D[j] / ( ( (Y[j] - 1) - D[j] ) * ( (Y[j] - 1) - D[j] ) * Y[j] * Y[j] );
					  					    }
					    				    else {
					    				      if (obsT[k] == time[j]){
					    					      c_kl_int = (Y[j] - D[j]) * (Y[j] - 1) * D[j] / ( ( (Y[j] - 1) - (D[j] - status[k]) ) * ( (Y[j] - 1) - (D[j] - status[k]) ) * Y[j] * Y[j] );
					    					    }
					    				      else { // obsT[k] < time[j]
					    					      c_kl_int = D[j] / ( (Y[j] - D[j]) * Y[j]); // (Y[j] - D[j]) * Y[j] * D[j] / ( (Y[j] - D[j]) * (Y[j] - D[j]) * Y[j] * Y[j] ) reduces to 
					    					    }
					        				}
												}
											else 		
												{//Reduced sample KM
													if (obsT[k] > time[j]){
					          				c_kl_int = D[j] / ( (Y[j] - 1 - D[j]) * (Y[j] - 1) );
					        				}
					        				else {
					          				if (obsT[k] == time[j]){
					            				c_kl_int = (D[j] - status[k]) / ( ( Y[j] - 1 - (D[j] - status[k]) ) * (Y[j] - 1) );
					          				}
					          				else { // obsT[k] < time[j]
					            				c_kl_int = D[j] / ( (Y[j] - D[j]) * Y[j]);
					          				}
					        				}
												}
											covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] = covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] + c_kl_int;
											// v2 = v2 + Sk[k] * Sk[k] * covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] / Sp[t*(*N) + k] / Sp[t*(*N) + k];
											// l > k case
											for (l=k+1;l<*N;l++){ // beginning of l-loop
												if (k==0){
													if (obsT[l]>time[j])
														{
															// decrease the number at risk because individual l was at risk at time[j] 
															na_k = D[j] / (Y[j] - 1);
															c_k_int = ( (Y[j] - 1) * D[j] ) / ( Y[j] * Y[j] * ( Y[j] - 1 - D[j] ) );
														}
													else{
															if (obsT[l]==time[j])
																{
																	// decrease the number of events if l was an event, and decrease the number at risk because l was in the risk set at time[j] 
																	na_k = (D[j] - status[l]) / (Y[j] - 1);
																	c_k_int = ( (Y[j] - 1) * D[j] ) / ( Y[j] * Y[j] * ( Y[j] - 1 - ( D[j] - status[l] ) ) );  
																}
															else{ // obsT[l]<time[j]
																	// do nothing since l is no longer in the risk set 
																	na_k = D[j] / Y[j];
																	c_k_int = D[j] / ( Y[j] * ( Y[j] - D[j] ) ); // ( Y[j] * D[j]) / ( Y[j] * Y[j] * ( Y[j] - D[j] ) ) reduces to
																}
														}
													Sk[l] = Sk[l] * (1 - na_k);
													covSSk[l] = covSSk[l] + c_k_int;
													}
												// k!= 0
												if (obsT[k] > time[j]){
											    if (obsT[l] > time[j]){
									      		c_kl_int = ( (Y[j] - D[j]) * (Y[j] - 2) * D[j] ) / ( (Y[j] - 1 - D[j]) * (Y[j] - 1 - D[j]) * Y[j] * Y[j] );
									    		}
									    		else {
									      		if (obsT[l] == time[j]){
									      		  c_kl_int = ( (Y[j] - D[j]) * (Y[j] - 2) * D[j] ) / ( (Y[j] - 1 - D[j]) * (Y[j] - 1 - (D[j] - status[l])) * Y[j] * Y[j] );
									      		}
									      		else { // obsT[l] < time[j]
									      		  c_kl_int = ( (Y[j] - 1) * D[j] ) / ( (Y[j] - 1 - D[j]) * Y[j] * Y[j] ); // ( (Y[j] - D[j]) * (Y[j] - 1) * D[j] ) / ( (Y[j] - 1 - D[j]) * (Y[j] - D[j]) * Y[j] * Y[j] ) reduces to 
									      		}
									    		}
									  		}
									  		else {
									    		if (obsT[k] == time[j]){
									      		if (obsT[l] > time[j]){
									        		c_kl_int = ( (Y[j] - D[j]) * (Y[j] - 2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k])) * (Y[j] - 1 - D[j]) * Y[j] * Y[j] );
									      		}
									      		else {
									        		if (obsT[l] == time[j]){
									          		c_kl_int = ( (Y[j] - D[j]) * (Y[j] - 2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k])) * (Y[j] - 1 - (D[j] - status[l])) * Y[j] * Y[j] );
									        		}
									        		else { // obsT[l] < time[j]
									          		c_kl_int = ( (Y[j] - 1) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k])) * Y[j] * Y[j] ); // ( (Y[j] - D[j]) * (Y[j] - 1) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k])) * (Y[j] - D[j]) * Y[j] * Y[j] ) reduces to 
									        		}
				 		                }	
											    }
											    else { // obsT[k] < time[j]
											      if (obsT[l] > time[j]){
											        c_kl_int = ( (Y[j] - 1) * D[j] ) / ( (Y[j] - 1 - D[j]) * Y[j] * Y[j] ); // ( (Y[j] - D[j]) * (Y[j] - 1) * D[j] ) / ( (Y[j] - D[j]) * (Y[j] - 1 - D[j]) * Y[j] * Y[j] ) reduces to 
											      }
											      else {
											        if (obsT[l] == time[j]){
											          c_kl_int = ( (Y[j] - 1) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[l])) * Y[j] * Y[j] ); // ( (Y[j] - D[j]) * (Y[j] - 1) * D[j] ) / ( (Y[j] - D[j]) * (Y[j] - 1 - (D[j] - status[l])) * Y[j] * Y[j] ) reduces to 
											        }
											        else { // obsT[l] < time[j]
											          c_kl_int = D[j] / ( (Y[j] - D[j]) * Y[j] ); // ( (Y[j] - D[j]) * Y[j] * D[j] ) / ( (Y[j] - D[j]) * (Y[j] - D[j]) * Y[j] * Y[j] ) reduces to 
											        }
											      }
											    }
											  }
												covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] = covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] + c_kl_int;
												// v2 = v2 + 2 * Sk[k] * Sk[l] * covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] / Sp[t*(*N) + k] / Sp[t*(*N) + l];
												} // end of l-loop
											// Compute the estimates 
											// S_Oi[k] = (*N) * S - (*N-1) * Sk[k];
											// S_E[t] = S_E[t] + S_Oi[k]/Sp[t*(*N) + k];
											// v3 = v3 + S * Sk[k] * covSSk[k] / Sp[t*(*N) + k];
											} // end of k-loop
											j++;
										} // end of while loop over j
									// j==i_t
									S = S * (1 - D[j] / Y[j]);
									varS = varS + D[j] / ( (Y[j] - D[j]) * Y[j] );
									for (k=0;k<*N;k++){
										v1 = v1 + 1 / Sp[t*(*N) + k];
										if (k==0) {
											if (obsT[k]>time[i_t])
												{
													// decrease the number at risk because individual k was at risk at time[i_t]
													na_k = D[i_t] / (Y[i_t] - 1);
													c_k_int = ( (Y[i_t] - 1) * D[i_t] ) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - 1 - D[i_t] ) );
												}
											else{
													if (obsT[k]==time[i_t])
														{
															// decrease the number of events if k was an event, and decrease the number at risk because k was in the risk set at time[i_t]
															na_k = (D[i_t] - status[k]) / (Y[i_t] - 1);
															c_k_int = ( (Y[i_t] - 1) * D[i_t] ) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - 1 - ( D[i_t] - status[k] ) ) );  
														}
													else{ // obsT[k]<time[i_t]
															// do nothing since k is no longer in the risk set 
															na_k = D[i_t] / Y[i_t];
															c_k_int = D[i_t] / ( Y[i_t] * ( Y[i_t] - D[i_t] ) ); // ( Y[i_t] * D[i_t]) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - D[i_t] ) ) reduces to
														}
												}
											Sk[k] = Sk[k] * (1 - na_k);
											covSSk[k] = covSSk[k] + c_k_int;
											} // end of k==0 case
										// l == k case
										if (*varSkMethod == 1)
											{//precise
												if (obsT[k] > time[i_t]){
					 	              c_kl_int = (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] / ( ( (Y[i_t] - 1) - D[i_t] ) * ( (Y[i_t] - 1) - D[i_t] ) * Y[i_t] * Y[i_t] );
		      					    }
		        				    else {
		        				      if (obsT[k] == time[i_t]){
		        					      c_kl_int = (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] / ( ( (Y[i_t] - 1) - (D[i_t] - status[k]) ) * ( (Y[i_t] - 1) - (D[i_t] - status[k]) ) * Y[i_t] * Y[i_t] );
		        					    }
		        				      else { // obsT[k] < time[i_t]
		        					      c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t]) * Y[i_t]); // (Y[i_t] - D[i_t]) * Y[i_t] * D[i_t] / ( (Y[i_t] - D[i_t]) * (Y[i_t] - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to 
		        					    }
		            				}
											}
										else 		
											{//Reduced sample KM
												if (obsT[k] > time[i_t]){
		              				c_kl_int = D[i_t] / ( (Y[i_t] - 1 - D[i_t]) * (Y[i_t] - 1) );
		            				}
		            				else {
		              				if (obsT[k] == time[i_t]){
		                				c_kl_int = (D[i_t] - status[k]) / ( ( Y[i_t] - 1 - (D[i_t] - status[k]) ) * (Y[i_t] - 1) );
		              				}
		              				else { // obsT[k] < time[i_t]
		                				c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t]) * Y[i_t]);
		              				}
		            				}
											}
										covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] = covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] + c_kl_int;
										v2 = v2 + Sk[k] * Sk[k] * covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] / Sp[t*(*N) + k] / Sp[t*(*N) + k];
										// l > k case
										for (l=k+1;l<*N;l++){ // beginning of l-loop
											if (k==0){
												if (obsT[l]>time[i_t])
													{
														// decrease the number at risk because individual l was at risk at time[i_t] 
														na_k = D[i_t] / (Y[i_t] - 1);
														c_k_int = ( (Y[i_t] - 1) * D[i_t] ) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - 1 - D[i_t] ) );
													}
												else{
														if (obsT[l]==time[i_t])
															{
																// decrease the number of events if l was an event, and decrease the number at risk because l was in the risk set at time[i_t] 
																na_k = (D[i_t] - status[l]) / (Y[i_t] - 1);
																c_k_int = ( (Y[i_t] - 1) * D[i_t] ) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - 1 - ( D[i_t] - status[l] ) ) );  
															}
														else{ // obsT[l]<time[i_t]
																// do nothing since l is no longer in the risk set 
																na_k = D[i_t] / Y[i_t];
																c_k_int = D[i_t] / ( Y[i_t] * ( Y[i_t] - D[i_t] ) ); // ( Y[i_t] * D[i_t]) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - D[i_t] ) ) reduces to
															}
													}
												Sk[l] = Sk[l] * (1 - na_k);
												covSSk[l] = covSSk[l] + c_k_int;
												}
											// k!= 0
											if (obsT[k] > time[i_t]){
								        if (obsT[l] > time[i_t]){
						          		c_kl_int = ( (Y[i_t] - D[i_t]) * (Y[i_t] - 2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] );
						        		}
						        		else {
						          		if (obsT[l] == time[i_t]){
						          		  c_kl_int = ( (Y[i_t] - D[i_t]) * (Y[i_t] - 2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * (Y[i_t] - 1 - (D[i_t] - status[l])) * Y[i_t] * Y[i_t] );
						          		}
						          		else { // obsT[l] < time[i_t]
						          		  c_kl_int = ( (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * (Y[i_t] - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to 
						          		}
						        		}
						      		}
						      		else {
						        		if (obsT[k] == time[i_t]){
						          		if (obsT[l] > time[i_t]){
						            		c_kl_int = ( (Y[i_t] - D[i_t]) * (Y[i_t] - 2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k])) * (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] );
						          		}
						          		else {
						            		if (obsT[l] == time[i_t]){
						              		c_kl_int = ( (Y[i_t] - D[i_t]) * (Y[i_t] - 2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k])) * (Y[i_t] - 1 - (D[i_t] - status[l])) * Y[i_t] * Y[i_t] );
						            		}
						            		else { // obsT[l] < time[i_t]
						              		c_kl_int = ( (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k])) * Y[i_t] * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k])) * (Y[i_t] - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to
						            		}
			 		                }	
								        }
								        else { // obsT[k] < time[i_t]
								          if (obsT[l] > time[i_t]){
								            c_kl_int = ( (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to 
								          }
								          else {
								            if (obsT[l] == time[i_t]){
								              c_kl_int = ( (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[l])) * Y[i_t] * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1 - (D[i_t] - status[l])) * Y[i_t] * Y[i_t] ) reduces to
								            }
								            else { // obsT[l] < time[i_t]
								              c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t]) * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * Y[i_t] * D[i_t] ) / ( (Y[i_t] - D[i_t]) * (Y[i_t] - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to 
								            }
								          }
								        }
								      }
											covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] = covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] + c_kl_int;
											v2 = v2 + 2 * Sk[k] * Sk[l] * covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] / Sp[t*(*N) + k] / Sp[t*(*N) + l];
											} // end of l-loop
										// Compute the estimates 
										S_Oi[k] = (*N) * S - (*N-1) * Sk[k];
										S_E[t] = S_E[t] + S_Oi[k]/Sp[t*(*N) + k];
										v3 = v3 + S * Sk[k] * covSSk[k] / Sp[t*(*N) + k];
										} // end of k-loop
								}						
							else { // at event times
									S = S * (1 - D[i_t] / Y[i_t]);
									varS = varS + D[i_t] / ( (Y[i_t] - D[i_t]) * Y[i_t] );
									for (k=0;k<*N;k++){
										v1 = v1 + 1 / Sp[t*(*N) + k];
										if (k==0) {
											if (obsT[k]>time[i_t])
												{
													// decrease the number at risk because individual k was at risk at time[i_t] 
													na_k = D[i_t] / (Y[i_t] - 1);
													c_k_int = ( (Y[i_t] - 1) * D[i_t] ) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - 1 - D[i_t] ) );
												}
											else{
													if (obsT[k]==time[i_t])
														{
															// decrease the number of events if k was an event, and decrease the number at risk because k was in the risk set at time[i_t] 
															na_k = (D[i_t] - status[k]) / (Y[i_t] - 1);
															c_k_int = ( (Y[i_t] - 1) * D[i_t] ) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - 1 - ( D[i_t] - status[k] ) ) );  
														}
													else{ // obsT[k]<time[i_t]
															// do nothing since k is no longer in the risk set 
															na_k = D[i_t] / Y[i_t];
															c_k_int = D[i_t] / ( Y[i_t] * ( Y[i_t] - D[i_t] ) ); // ( Y[i_t] * D[i_t]) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - D[i_t] ) ) reduces to
														}
												}
											Sk[k] = Sk[k] * (1 - na_k);
											covSSk[k] = covSSk[k] + c_k_int;
											} // end of k==0 case
										// l == k case
										if (*varSkMethod == 1)
											{//precise
												if (obsT[k] > time[i_t]){
					 	              c_kl_int = (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] / ( ( (Y[i_t] - 1) - D[i_t] ) * ( (Y[i_t] - 1) - D[i_t] ) * Y[i_t] * Y[i_t] );
		      					    }
		        				    else {
		        				      if (obsT[k] == time[i_t]){
		        					      c_kl_int = (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] / ( ( (Y[i_t] - 1) - (D[i_t] - status[k]) ) * ( (Y[i_t] - 1) - (D[i_t] - status[k]) ) * Y[i_t] * Y[i_t] );
		        					    }
		        				      else { // obsT[k] < time[i_t]
		        					      c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t]) * Y[i_t]); // (Y[i_t] - D[i_t]) * Y[i_t] * D[i_t] / ( (Y[i_t] - D[i_t]) * (Y[i_t] - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to
		        					    }
		            				}
											}
										else 		
											{//Reduced sample KM
												if (obsT[k] > time[i_t]){
		              				c_kl_int = D[i_t] / ( (Y[i_t] - 1 - D[i_t]) * (Y[i_t] - 1) );
		            				}
		            				else {
		              				if (obsT[k] == time[i_t]){
		                				c_kl_int = (D[i_t] - status[k]) / ( ( Y[i_t] - 1 - (D[i_t] - status[k]) ) * (Y[i_t] - 1) );
		              				}
		              				else { // obsT[k] < time[i_t]
		                				c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t]) * Y[i_t]);
		              				}
		            				}
											}
										covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] = covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] + c_kl_int;
										v2 = v2 + Sk[k] * Sk[k] * covSkSl[(*N)*(*N + 1)/2 - (*N-k)*(*N+1-k)/2] / Sp[t*(*N) + k] / Sp[t*(*N) + k];
										// l > k case
										for (l=k+1;l<*N;l++){ // beginning of l-loop
											if (k==0){
												if (obsT[l]>time[i_t])
													{
														// decrease the number at risk because individual l was at risk at time[i_t] 
														na_k = D[i_t] / (Y[i_t] - 1);
														c_k_int = ( (Y[i_t] - 1) * D[i_t] ) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - 1 - D[i_t] ) );
													}
												else{
														if (obsT[l]==time[i_t])
															{
																// decrease the number of events if l was an event, and decrease the number at risk because l was in the risk set at time[i_t] 
																na_k = (D[i_t] - status[l]) / (Y[i_t] - 1);
																c_k_int = ( (Y[i_t] - 1) * D[i_t] ) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - 1 - ( D[i_t] - status[l] ) ) );  
															}
														else{ // obsT[l]<time[i_t]
																// do nothing since l is no longer in the risk set 
																na_k = D[i_t] / Y[i_t];
																c_k_int = D[i_t] / ( Y[i_t] * ( Y[i_t] - D[i_t] ) ); // ( Y[i_t] * D[i_t]) / ( Y[i_t] * Y[i_t] * ( Y[i_t] - D[i_t] ) ) reduces to
															}
													}
												Sk[l] = Sk[l] * (1 - na_k);
												covSSk[l] = covSSk[l] + c_k_int;
												}
											// k!= 0
											if (obsT[k] > time[i_t]){
								        if (obsT[l] > time[i_t]){
						          		c_kl_int = ( (Y[i_t] - D[i_t]) * (Y[i_t] - 2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] );
						        		}
						        		else {
						          		if (obsT[l] == time[i_t]){
						          		  c_kl_int = ( (Y[i_t] - D[i_t]) * (Y[i_t] - 2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * (Y[i_t] - 1 - (D[i_t] - status[l])) * Y[i_t] * Y[i_t] );
						          		}
						          		else { // obsT[l] < time[i_t]
						          		  c_kl_int = ( (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * (Y[i_t] - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to 
						          		}
						        		}
						      		}
						      		else {
						        		if (obsT[k] == time[i_t]){
						          		if (obsT[l] > time[i_t]){
						            		c_kl_int = ( (Y[i_t] - D[i_t]) * (Y[i_t] - 2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k])) * (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] );
						          		}
						          		else {
						            		if (obsT[l] == time[i_t]){
						              		c_kl_int = ( (Y[i_t] - D[i_t]) * (Y[i_t] - 2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k])) * (Y[i_t] - 1 - (D[i_t] - status[l])) * Y[i_t] * Y[i_t] );
						            		}
						            		else { // obsT[l] < time[i_t]
						              		c_kl_int = ( (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k])) * Y[i_t] * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k])) * (Y[i_t] - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to
						            		}
			 		                }	
								        }
								        else { // obsT[k] < time[i_t]
								          if (obsT[l] > time[i_t]){
								            c_kl_int = ( (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1 - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to 
								          }
								          else {
								            if (obsT[l] == time[i_t]){
								              c_kl_int = ( (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[l])) * Y[i_t] * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1) * D[i_t] ) / ( (Y[i_t] - D[i_t]) * (Y[i_t] - 1 - (D[i_t] - status[l])) * Y[i_t] * Y[i_t] ) reduces to
								            }
								            else { // obsT[l] < time[i_t]
								              c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t]) * Y[i_t] ); // ( (Y[i_t] - D[i_t]) * Y[i_t] * D[i_t] ) / ( (Y[i_t] - D[i_t]) * (Y[i_t] - D[i_t]) * Y[i_t] * Y[i_t] ) reduces to
								            }
								          }
								        }
								      }
											covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] = covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] + c_kl_int;
											v2 = v2 + 2 * Sk[k] * Sk[l] * covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)] / Sp[t*(*N) + k] / Sp[t*(*N) + l];
											} // end of l-loop
										// Compute the estimates 
										S_Oi[k] = (*N) * S - (*N-1) * Sk[k];
										S_E[t] = S_E[t] + S_Oi[k]/Sp[t*(*N) + k];
										v3 = v3 + S * Sk[k] * covSSk[k] / Sp[t*(*N) + k];
										} // end of k-loop
								}
						}
				}
			//printf("\n t \n");printf("%d ",t);printf("\n i_t "); printf("%d ", i_t);
			//printf("\n S. ");printf("%f ",S);
			//printf("\n time je ");printf("%f ", time[t]);			
			/*for (k=0;k<*N;k++){
				printf("\n k = ");printf("%d ",k);
				printf("\n Sk[k] = ");printf("%f ",Sk[k]);
				printf("\n S_Oi[k] = ");printf("%f \n",S_Oi[k]);
				for (l=k;l<*N;l++){
					//printf(" l = ");printf("%d ",l);
					printf("%f ",covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)]);
					}		
			//printf("\n S_Oi[k] ");printf("%f ",S_Oi[k]);
			//printf("\n Sp[k] ");printf("%f ",Sp[t*(*N) + k]);
			//printf("\n S_Oi[k] / Sp[k] ");printf("%f ",S_Oi[k] / Sp[t*(*N) + k]);
			}*/
			S_E[t] = S_E[t] / (*N);
			//printf("\n varS je ");printf("%f ", S * S * varS);
			//printf("\n t je ");printf("%d ", t);
			/*for (k=0;k<*N;k++){
			//printf("\n Sk ");printf("%f ",Sk[k]);
			printf("\n covSSk ");printf("%f ",covSSk[k]);
			}*/
			v2 = v2 * (*N-1) * (*N-1) / (*N) / (*N);
			v3 = - v3 * v1 * 2 * (*N-1) / (*N);
			v1 = v1 * v1 * S * S * varS;
			//printf("\n v1 ");printf("%f ", v1);
			//printf("\n v2 ");printf("%f ", v2);
			//printf("\n v3 ");printf("%f ", -v3);
			varS_E[t] = v1 + v2 + v3;
		}
	}
/*
#include <stdio.h>
int main()
{	
	double Y[] = {10, 8, 7, 5, 3, 2};
	double D[] = {1, 1, 2, 1, 1, 1};
	double time[] = { 70, 458, 475, 861, 1196, 1360};
	double obsT[] = {70, 70, 458, 475, 475, 861, 980, 1196, 1360, 3652};
	double status[] = {1, 0, 1, 1, 1, 1, 0, 1, 1, 0};
	//int neval = 10;
	int N = 10;
	int NT = 6;
	//int ind[] = {-1, 0, 1, 2, 2, 3, 3, 4, 5, 5};
	//double evalT[] = {50, 70, 458, 475, 800, 861, 980, 1196, 1360, 3652};	
	//double Sp[] = {0.9965768, 0.9993580, 0.9996044, 0.9950994, 0.9973064, 0.9957560, 0.9992099, 0.9983667, 0.9834300, 0.9853559, 0.9952140, 0.9991014, 0.9994289, 0.9931458, 0.9962310, 0.9940634, 0.9988763, 0.9977141, 0.9768680, 0.9795585, 0.9688255, 0.9945278, 0.9957526, 0.9531374, 0.9749952, 0.9622891, 0.9924566, 0.9850181, 0.8584975, 0.8673584, 0.9674023, 0.9943339, 0.9956101, 0.9513747, 0.9740650, 0.9609866, 0.9921894, 0.9844661, 0.8530743, 0.8623909, 0.9407879, 0.9909505, 0.9927721, 0.9159922, 0.9548499, 0.9360654, 0.9870759, 0.9722998, 0.7485556, 0.7693454, 0.9365613, 0.9902314, 0.9922159, 0.9094313, 0.9503742, 0.9313205, 0.9863575, 0.9699216, 0.7292449, 0.7522086, 0.9283707, 0.9887126, 0.9911156, 0.8970927, 0.9417454, 0.9221331, 0.9849575, 0.9651351, 0.6925474, 0.7199840, 0.9132076, 0.9861517, 0.9891135, 0.8725600, 0.9260404, 0.9044571, 0.9815183, 0.9562944, 0.6306977, 0.6663558, 0.9005569, 0.9843874, 0.9877244, 0.8545227, 0.9132742, 0.8903809, 0.9787590, 0.9491318, 0.5870748, 0.6205567, 0.7174663, 0.9505583, 0.9595669, 0.4885878, 0.7114414, 0.6915124, 0.9272684, 0.8465829, 0.1478734, 0.1646298};
	// at jump times only
	//int neval = 6;
	//int ind[] = {0, 1, 2, 3, 4, 5};
	//double evalT[] = {70, 458, 475, 861, 1196, 1360};	
	//double Sp[]={0.9952140, 0.9991014, 0.9994289, 0.9931458, 0.9962310, 0.9940634, 0.9988763, 0.9977141, 0.9768680, 0.9795585, 0.9688255, 0.9945278, 0.9957526, 0.9531374, 0.9749952, 0.9622891, 0.9924566, 0.9850181, 0.8584975, 0.8673584, 0.9674023, 0.9943339, 0.9956101, 0.9513747, 0.9740650, 0.9609866, 0.9921894, 0.9844661, 0.8530743, 0.8623909, 0.9365613, 0.9902314, 0.9922159, 0.9094313, 0.9503742, 0.9313205, 0.9863575, 0.9699216, 0.7292449, 0.7522086, 0.9132076, 0.9861517, 0.9891135, 0.8725600, 0.9260404, 0.9044571, 0.9815183, 0.9562944, 0.6306977, 0.6663558, 0.9005569, 0.9843874, 0.9877244, 0.8545227, 0.9132742, 0.8903809, 0.9787590, 0.9491318, 0.5870748, 0.6205567};	
	// special cases
	int neval = 1;
	int ind[] = {2};
	double evalT[] = {475};
	double Sp[] = {0.9674023, 0.9943339, 0.9956101, 0.9513747, 0.9740650, 0.9609866, 0.9921894, 0.9844661, 0.8530743, 0.8623909};

	double S_E[neval];
	double varS_E[neval];	
	//int varSkMethod = 1; //"precise"
	int varSkMethod = 2; //"reduced sample (KM)"
	double covSkSl[N*(N+1)/2];

	double *Y_p = &Y[0];
	double *D_p = &D[0]; 
	double *time_p = &time[0];	
	double *obsT_p = &obsT[0];
	double *status_p = &status[0]; 
	int *neval_p;
	neval_p = &neval;
	int *N_p;
	N_p = &N;
	int *NT_p;
	NT_p = &NT;
	int *ind_p = &ind[0];
	double *Sp_p = &Sp[0];
	double *S_E_p = &S_E[0];
	double *varS_E_p = &varS_E[0];
	int *varSkMethod_p;
	varSkMethod_p = &varSkMethod;	
	double *covSkSl_p = &covSkSl[0];

	int i;

	//pseudoWH_pre2(Y_p, D_p, time_p, obsT_p, status_p, S_E_p, varS_E_p, Sp_p, N_p, NT_p, neval_p, ind_p, varSkMethod_p, covSkSl_p);
	
	
	printf("\nPo evaluaciji.\n");	
	printf("Vrednost S_E je:\n");
	for (i=0;i<neval;i++){
		printf("%d ",i); 
		printf("%f ",S_E[i]);
	}
	printf("\nVrednost varS_E je:\n");
	for(i=0;i<neval;i++){
		printf("%f ",varS_E[i]);		
		//printf("%f\n", *(varS_p+i));
	}
	
	return 0;
}
*/
