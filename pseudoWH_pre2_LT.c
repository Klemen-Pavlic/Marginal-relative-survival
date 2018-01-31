/* 
	Author: Klemen Pavlic 
*/

void pseudoWH_pre2_LT(double *Y, double *D, double *C, double *time, double *obsT, double *status, double *S_E, double *varS_E, double *Sp, int *N, int *NT, int *neval, int *ind, int *varSkMethod, double *covSkSl){
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
									S = S * (1 - D[j] / (Y[j] - C[j]/2));
									varS = varS + D[j] / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) );
									for (k=0;k<*N;k++){
										//v1 = v1 + 1 / Sp[t*(*N) + k];
										if (k==0) {
											if (obsT[k]>time[j])
												{
													// decrease the number at risk because individual k was at risk at time[j] 
													na_k = D[j] / (Y[j] - 1 - C[j]/2);
													c_k_int = ( (Y[j] - 1 - C[j]/2) * D[j] ) / ( (Y[j] - C[j]/2) * (Y[j] - C[j]/2) * ( Y[j] - 1 - D[j] - C[j]/2) );
												}
											else{
													if (obsT[k]==time[j])
														{
															// decrease the number of events if k was an event, and decrease the number at risk because k was in the risk set at time[j]
															na_k = (D[j] - status[k]) / (Y[j] - 1 - (C[j] - (1-status[k]))/2);
															c_k_int = ( (Y[j] - 1 - (C[j] - (1-status[k]))/2) * D[j] ) / ( (Y[j] - C[j]/2) * (Y[j] - C[j]/2) * ( Y[j] - 1 - ( D[j] - status[k] ) - (C[j] - (1-status[k]))/2 ) );  
														}
													else{ // obsT[k]<time[i_t]
															// do nothing since k is no longer in the risk set 
															na_k = D[j] / (Y[j] - C[j]/2);
															c_k_int = D[j] / ( (Y[j] - C[j]/2) * ( Y[j] - D[j] - C[j]/2 ) ); // ( (Y[j] - C[j]/2) * D[j]) / ( (Y[j] - C[j]/2) * (Y[j] - C[j]/2) * ( Y[j] - D[j] - C[j]/2 ) ) reduces to
														}
												}
											Sk[k] = Sk[k] * (1 - na_k);
											covSSk[k] = covSSk[k] + c_k_int;
											} // end of k==0 case
										// l == k case
										if (*varSkMethod == 1)
											{//precise
												if (obsT[k] > time[j]){
					 	              c_kl_int = (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - C[j]/2) * D[j] / ( ( (Y[j] - 1) - D[j] - C[j]/2 ) * ( (Y[j] - 1) - D[j] - C[j]/2 ) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) );
		      					    }
		        				    else {
		        				      if (obsT[k] == time[j]){
		        					      c_kl_int = (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - (C[j] - (1- status[k]))/2) * D[j] / ( ( (Y[j] - 1) - (D[j] - status[k]) - (C[j] - (1-status[k]))/2 ) * ( (Y[j] - 1) - (D[j] - status[k]) - (C[j] - (1-status[k]))/2 ) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) );
		        					    }
		        				      else { // obsT[k] < time[j]
		        					      c_kl_int = D[j] / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2)); // (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) * D[j] / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ) reduces to 
		        					    }
		            				}
											}
										else 		
											{//Reduced sample KM
												if (obsT[k] > time[j]){
		              				c_kl_int = D[j] / ( (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - 1 - C[j]/2) );
		            				}
		            				else {
		              				if (obsT[k] == time[j]){
		                				c_kl_int = (D[j] - status[k]) / ( ( Y[j] - 1 - (D[j] - status[k]) - (C[j] - (1-status[k]))/2 ) * (Y[j] - 1 - (C[j] - (1-status[k]))/2) );
		              				}
		              				else { // obsT[k] < time[j]
		                				c_kl_int = D[j] / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) );
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
														na_k = D[j] / (Y[j] - 1 - C[j]/2);
														c_k_int = ( (Y[j] - 1 - C[j]/2) * D[j] ) / ( (Y[j] - C[j]/2) * (Y[j] - C[j]/2) * ( Y[j] - 1 - D[j] - C[j]/2 ) );
													}
												else{
														if (obsT[l]==time[j])
															{
																// decrease the number of events if l was an event, and decrease the number at risk because l was in the risk set at time[j] 
																na_k = (D[j] - status[l]) / (Y[j] - 1 - (C[j] - (1-status[l]))/2);
																c_k_int = ( (Y[j] - 1 - (C[j] - (1-status[l]))/2) * D[j] ) / ( (Y[j] - C[j]/2) * (Y[j] - C[j]/2) * ( Y[j] - 1 - ( D[j] - status[l] ) - (C[j] - (1-status[l]))/2 ) );  
															}
														else{ // obsT[l]<time[j]
																// do nothing since l is no longer in the risk set 
																na_k = D[j] / (Y[j] - C[j]/2);
																c_k_int = D[j] / ( (Y[j] - C[j]/2) * ( Y[j] - D[j] - C[j]/2 ) ); // ( (Y[j] - C[j]/2) * D[j]) / ( (Y[j] - C[j]/2) * (Y[j] - C[j]/2) * ( Y[j] - D[j] - C[j]/2) ) reduces to
															}
													}
												Sk[l] = Sk[l] * (1 - na_k);
												covSSk[l] = covSSk[l] + c_k_int;
												}
											// k!= 0
											if (obsT[k] > time[j]){
								        if (obsT[l] > time[j]){
						          		c_kl_int = ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 2 - C[j]/2) * D[j] ) / ( (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) );
						        		}
						        		else {
						          		if (obsT[l] == time[j]){
						          		  c_kl_int = ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 2 - (C[j] - (1-status[l]))/2) * D[j] ) / ( (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - 1 - (D[j] - status[l]) - (C[j] - (1-status[l]))/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) );
						          		}
						          		else { // obsT[l] < time[j]
						          		  c_kl_int = ( (Y[j] - 1 - C[j]/2) * D[j] ) / ( (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ); // ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - C[j]/2) * D[j] ) / ( (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ) reduces to 
						          		}
						        		}
						      		}
						      		else {
						        		if (obsT[k] == time[j]){
						          		if (obsT[l] > time[j]){
						            		c_kl_int = ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 2 - (C[j] - (1-status[k]))/2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k]) - (C[j] - (1-status[k]))/2) * (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) );
						          		}
						          		else {
						            		if (obsT[l] == time[j]){
						              		c_kl_int = ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 2 - (C[j] - (1-status[k]) - (1-status[l]))/2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k]) - (C[j] - (1-status[k]))/2) * (Y[j] - 1 - (D[j] - status[l]) - (C[j] - (1-status[l]))/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) );
						            		}
						            		else { // obsT[l] < time[j]
						              		c_kl_int = ( (Y[j] - 1 - (C[j] - (1-status[k]))/2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k]) - (C[j] - (1-status[k]))/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ); // ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - (C[j] - (1-status[k]))/2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k]) - (C[j] - (1-status[k]))/2) * (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ) reduces to 
						            		}
			 		                }	
								        }
								        else { // obsT[k] < time[j]
								          if (obsT[l] > time[j]){
								            c_kl_int = ( (Y[j] - 1 - C[j]/2) * D[j] ) / ( (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ); // ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - C[j]/2) * D[j] ) / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ) reduces to 
								          }
								          else {
								            if (obsT[l] == time[j]){
								              c_kl_int = ( (Y[j] - 1 - (C[j] - (1-status[l]))/2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[l]) - (C[j] - (1-status[l]))/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ); // ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - (C[j] - (1-status[l]))/2) * D[j] ) / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - (D[j] - status[l]) - (C[j] - (1-status[l]))/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ) reduces to 
								            }
								            else { // obsT[l] < time[j]
								              c_kl_int = D[j] / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) ); // ( (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) * D[j] ) / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ) reduces to 
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
							S = S * (1 - D[i_t] / (Y[i_t] - C[i_t]/2));
							varS = varS + D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
							for (k=0;k<*N;k++){
								v1 = v1 + 1 / Sp[t*(*N) + k];
								if (k==0) {
									if (obsT[k]>time[i_t])
										{
											// decrease the number at risk because individual k was at risk at time[i_t]
											na_k = D[i_t] / (Y[i_t] - 1 - C[i_t]/2);
											c_k_int = ( (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - 1 - D[i_t] - C[i_t]/2 ) );
										}
									else{
											if (obsT[k]==time[i_t])
												{
													// decrease the number of events if k was an event, and decrease the number at risk because k was in the risk set at time[i_t] 
													na_k = (D[i_t] - status[k]) / (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2);
													c_k_int = ( (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) * D[i_t] ) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - 1 - ( D[i_t] - status[k] ) - (C[i_t] - (1-status[k]))/2 ) );  
												}
											else{ // obsT[k]<time[i_t]
													// do nothing since k is no longer in the risk set 
													na_k = D[i_t] / (Y[i_t] - C[i_t]/2);
													c_k_int = D[i_t] / ( (Y[i_t] - C[i_t]/2) * ( Y[i_t] - D[i_t] - C[i_t]/2 ) ); // ( (Y[i_t] - C[i_t]/2) * D[i_t]) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - D[i_t] - C[i_t]/2 ) ) reduces to
												}
										}
									Sk[k] = Sk[k] * (1 - na_k);
									covSSk[k] = covSSk[k] + c_k_int;
									} // end of k==0 case
								// l == k case
								if (*varSkMethod == 1)
									{//precise
										if (obsT[k] > time[i_t]){
			 	              c_kl_int = (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) * D[i_t] / ( ( (Y[i_t] - 1) - D[i_t] - C[i_t]/2 ) * ( (Y[i_t] - 1) - D[i_t] - C[i_t]/2 ) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
      					    }
        				    else {
        				      if (obsT[k] == time[i_t]){
        					      c_kl_int = (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) * D[i_t] / ( ( (Y[i_t] - 1) - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2 ) * ( (Y[i_t] - 1) - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2 ) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
        					    }
        				      else { // obsT[k] < time[i_t]
        					      c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to 
        					    }
            				}
									}
								else 		
									{//Reduced sample KM
										if (obsT[k] > time[i_t]){
              				c_kl_int = D[i_t] / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) );
            				}
            				else {
              				if (obsT[k] == time[i_t]){
                				c_kl_int = (D[i_t] - status[k]) / ( ( Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2 ) * (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) );
              				}
              				else { // obsT[k] < time[i_t]
                				c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2));
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
												na_k = D[i_t] / (Y[i_t] - 1 - C[i_t]/2);
												c_k_int = ( (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - 1 - D[i_t] - C[i_t]/2 ) );
											}
										else{
												if (obsT[l]==time[i_t])
													{
														// decrease the number of events if l was an event, and decrease the number at risk because l was in the risk set at time[i_t] 
														na_k = (D[i_t] - status[l]) / (Y[i_t] - 1 - (C[i_t] - (1-status[l]))/2);
														c_k_int = ( (Y[i_t] - 1 - (C[i_t] - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - 1 - ( D[i_t] - status[l] ) - (C[i_t] - (1-status[l]))/2 ) );  
													}
												else{ // obsT[l]<time[i_t]
														// do nothing since l is no longer in the risk set 
														na_k = D[i_t] / (Y[i_t] - C[i_t]/2);
														c_k_int = D[i_t] / ( (Y[i_t] - C[i_t]/2) * ( Y[i_t] - D[i_t] - C[i_t]/2 ) ); // ( (Y[i_t] - C[i_t]/2) * D[i_t]) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - D[i_t] - C[i_t]/2 ) ) reduces to
													}
											}
										Sk[l] = Sk[l] * (1 - na_k);
										covSSk[l] = covSSk[l] + c_k_int;
										}
									// k!= 0
									if (obsT[k] > time[i_t]){
						        if (obsT[l] > time[i_t]){
				          		c_kl_int = ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 2 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
				        		}
				        		else {
				          		if (obsT[l] == time[i_t]){
				          		  c_kl_int = ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 2 - (C[i_t] - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (D[i_t] - status[l]) - (C[i_t] - (1-status[l]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
				          		}
				          		else { // obsT[l] < time[i_t]
				          		  c_kl_int = ( (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to 
				          		}
				        		}
				      		}
				      		else {
				        		if (obsT[k] == time[i_t]){
				          		if (obsT[l] > time[i_t]){
				            		c_kl_int = ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 2 - (C[i_t] - (1-status[k]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2) * (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
				          		}
				          		else {
				            		if (obsT[l] == time[i_t]){
				              		c_kl_int = ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 2 - (C[i_t] - (1-status[k]) - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2) * (Y[i_t] - 1 - (D[i_t] - status[l]) - (C[i_t] - (1-status[l]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
				            		}
				            		else { // obsT[l] < time[i_t]
				              		c_kl_int = ( (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2) * (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to 
				            		}
	 		                }	
						        }
						        else { // obsT[k] < time[i_t]
						          if (obsT[l] > time[i_t]){
						            c_kl_int = ( (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to 
						          }
						          else {
						            if (obsT[l] == time[i_t]){
						              c_kl_int = ( (Y[i_t] - 1 - (C[i_t] - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[l]) - (C[i_t] - (1-status[l]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (C[i_t] - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (D[i_t] - status[l]) - (C[i_t] - (1-status[l]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to 
						            }
						            else { // obsT[l] < time[i_t]
						              c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to
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
											S = S * (1 - D[j] / (Y[j] - C[j]/2));
											varS = varS + D[j] / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) );
											for (k=0;k<*N;k++){
											//v1 = v1 + 1 / Sp[t*(*N) + k];
											if (k==0) {
												if (obsT[k]>time[j])
													{
														// decrease the number at risk because individual k was at risk at time[j] 
														na_k = D[j] / (Y[j] - 1 - C[j]/2);
														c_k_int = ( (Y[j] - 1 - C[j]/2) * D[j] ) / ( (Y[j] - C[j]/2) * (Y[j] - C[j]/2) * ( Y[j] - 1 - D[j] - C[j]/2 ) );
													}
												else{
														if (obsT[k]==time[j])
															{
																// decrease the number of events if k was an event, and decrease the number at risk because k was in the risk set at time[j] 
																na_k = (D[j] - status[k]) / (Y[j] - 1 - (C[j] - (1-status[k]))/2);
																c_k_int = ( (Y[j] - 1 - (C[j] - (1-status[k]))/2) * D[j] ) / ( (Y[j] - C[j]/2) * (Y[j] - C[j]/2) * ( Y[j] - 1 - ( D[j] - status[k] ) - (C[j] - (1-status[k]))/2 ) );  
															}
														else{ // obsT[k]<time[j]
																// do nothing since k is no longer in the risk set 
																na_k = D[j] / (Y[j] - C[j]/2);
																c_k_int = D[j] / ( (Y[j] - C[j]/2) * ( Y[j] - D[j] - C[j]/2 ) ); // ( (Y[j] -C[j]/2) * D[j]) / ( (Y[j] - C[j]/2) * (Y[j] - C[j]/2) * ( Y[j] - D[j] - C[j]/2 ) ) reduces to
															}
													}
												Sk[k] = Sk[k] * (1 - na_k);
												covSSk[k] = covSSk[k] + c_k_int;
												} // end of k==0 case
											// l == k case
											if (*varSkMethod == 1)
												{//precise
													if (obsT[k] > time[j]){
						 	              c_kl_int = (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - C[j]/2) * D[j] / ( ( (Y[j] - 1) - D[j] - C[j]/2 ) * ( (Y[j] - 1) - D[j] - C[j]/2 ) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) );
					  					    }
					    				    else {
					    				      if (obsT[k] == time[j]){
					    					      c_kl_int = (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - (C[j] - (1-status[k]))/2) * D[j] / ( ( (Y[j] - 1) - (D[j] - status[k]) - (C[j] - (1-status[k]))/2 ) * ( (Y[j] - 1) - (D[j] - status[k]) - (C[j] - (1-status[k]))/2 ) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) );
					    					    }
					    				      else { // obsT[k] < time[j]
					    					      c_kl_int = D[j] / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) ); // (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) * D[j] / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ) reduces to 
					    					    }
					        				}
												}
											else 		
												{//Reduced sample KM
													if (obsT[k] > time[j]){
					          				c_kl_int = D[j] / ( (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - 1 - C[j]/2) );
					        				}
					        				else {
					          				if (obsT[k] == time[j]){
					            				c_kl_int = (D[j] - status[k]) / ( ( Y[j] - 1 - (D[j] - status[k]) - (C[j] - (1-status[k]))/2 ) * (Y[j] - 1 - (C[j] - (1-status[k]))/2) );
					          				}
					          				else { // obsT[k] < time[j]
					            				c_kl_int = D[j] / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) );
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
															na_k = D[j] / (Y[j] - 1 - C[j]/2);
															c_k_int = ( (Y[j] - 1 - C[j]/2) * D[j] ) / ( (Y[j] - C[j]/2) * (Y[j] - C[j]/2) * ( Y[j] - 1 - D[j] - C[j]/2 ) );
														}
													else{
															if (obsT[l]==time[j])
																{
																	// decrease the number of events if l was an event, and decrease the number at risk because l was in the risk set at time[j] 
																	na_k = (D[j] - status[l]) / (Y[j] - 1 - (C[j] - (1-status[l]))/2);
																	c_k_int = ( (Y[j] - 1 - (C[j] - (1-status[l]))/2) * D[j] ) / ( (Y[j] - C[j]/2) * (Y[j] - C[j]/2) * ( Y[j] - 1 - ( D[j] - status[l] ) - (C[j] - (1-status[l]))/2 ) );  
																}
															else{ // obsT[l]<time[j]
																	// do nothing since l is no longer in the risk set 
																	na_k = D[j] / (Y[j] - C[j]/2);
																	c_k_int = D[j] / ( (Y[j] - C[j]/2) * ( Y[j] - D[j] - C[j]/2 ) ); // ( (Y[j] - C[j]/2) * D[j]) / ( (Y[j] - C[j]/2) * (Y[j] - C[j]/2) * ( Y[j] - D[j] - C[j]/2 ) ) reduces to
																}
														}
													Sk[l] = Sk[l] * (1 - na_k);
													covSSk[l] = covSSk[l] + c_k_int;
													}
												// k!= 0
												if (obsT[k] > time[j]){
											    if (obsT[l] > time[j]){
									      		c_kl_int = ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 2 - C[j]/2) * D[j] ) / ( (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) );
									    		}
									    		else {
									      		if (obsT[l] == time[j]){
									      		  c_kl_int = ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 2 - (C[j] - (1-status[l]))/2) * D[j] ) / ( (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - 1 - (D[j] - status[l]) - (C[j] - (1-status[l]))/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) );
									      		}
									      		else { // obsT[l] < time[j]
									      		  c_kl_int = ( (Y[j] - 1 - C[j]/2) * D[j] ) / ( (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ); // ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - C[j]/2) * D[j] ) / ( (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ) reduces to 
									      		}
									    		}
									  		}
									  		else {
									    		if (obsT[k] == time[j]){
									      		if (obsT[l] > time[j]){
									        		c_kl_int = ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 2 - (C[j] - (1-status[k]))/2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k]) - (C[j] - (1-status[k]))/2) * (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) );
									      		}
									      		else {
									        		if (obsT[l] == time[j]){
									          		c_kl_int = ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 2 - (C[j] - (1-status[k]) - (1-status[l]))/2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k]) - (C[j] - (1-status[k]))/2) * (Y[j] - 1 - (D[j] - status[l]) - (C[j] - (1-status[l]))/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) );
									        		}
									        		else { // obsT[l] < time[j]
									          		c_kl_int = ( (Y[j] - 1 - (C[j] - (1-status[k]))/2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k]) - (C[j] - (1-status[k]))/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ); // ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - (C[j] - (1-status[k]))/2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[k]) - (C[j] - (1-status[k]))/2) * (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ) reduces to 
									        		}
				 		                }	
											    }
											    else { // obsT[k] < time[j]
											      if (obsT[l] > time[j]){
											        c_kl_int = ( (Y[j] - 1 - C[j]/2) * D[j] ) / ( (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ); // ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - C[j]/2) * D[j] ) / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ) reduces to 
											      }
											      else {
											        if (obsT[l] == time[j]){
											          c_kl_int = ( (Y[j] - 1 - (C[j] - (1-status[l]))/2) * D[j] ) / ( (Y[j] - 1 - (D[j] - status[l]) - (C[j] - (1-status[l]))/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ); // ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - (C[j] - (1-status[l]))/2) * D[j] ) / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - 1 - (D[j] - status[l]) - (C[j] - (1-status[l]))/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ) reduces to 
											        }
											        else { // obsT[l] < time[j]
											          c_kl_int = D[j] / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) ); // ( (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) * D[j] ) / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) * (Y[j] - C[j]/2) ) reduces to 
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
									S = S * (1 - D[j] / (Y[j] - C[j]/2));
									varS = varS + D[j] / ( (Y[j] - D[j] - C[j]/2) * (Y[j] - C[j]/2) );
									for (k=0;k<*N;k++){
										v1 = v1 + 1 / Sp[t*(*N) + k];
										if (k==0) {
											if (obsT[k]>time[i_t])
												{
													// decrease the number at risk because individual k was at risk at time[i_t]
													na_k = D[i_t] / (Y[i_t] - 1 - C[i_t]/2);
													c_k_int = ( (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - 1 - D[i_t] - C[i_t]/2 ) );
												}
											else{
													if (obsT[k]==time[i_t])
														{
															// decrease the number of events if k was an event, and decrease the number at risk because k was in the risk set at time[i_t]
															na_k = (D[i_t] - status[k]) / (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2);
															c_k_int = ( (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) * D[i_t] ) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - 1 - ( D[i_t] - status[k] ) - (C[i_t] - (1-status[k]))/2 ) );  
														}
													else{ // obsT[k]<time[i_t]
															// do nothing since k is no longer in the risk set 
															na_k = D[i_t] / (Y[i_t] - C[i_t]/2);
															c_k_int = D[i_t] / ( (Y[i_t] - C[i_t]/2) * ( Y[i_t] - D[i_t] - C[i_t]/2 ) ); // ( (Y[i_t] - C[i_t]/2) * D[i_t]) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - D[i_t] - C[i_t]/2 ) ) reduces to
														}
												}
											Sk[k] = Sk[k] * (1 - na_k);
											covSSk[k] = covSSk[k] + c_k_int;
											} // end of k==0 case
										// l == k case
										if (*varSkMethod == 1)
											{//precise
												if (obsT[k] > time[i_t]){
					 	              c_kl_int = (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) * D[i_t] / ( ( (Y[i_t] - 1) - D[i_t] - C[i_t]/2 ) * ( (Y[i_t] - 1) - D[i_t] - C[i_t]/2 ) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
		      					    }
		        				    else {
		        				      if (obsT[k] == time[i_t]){
		        					      c_kl_int = (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) * D[i_t] / ( ( (Y[i_t] - 1) - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2 ) * ( (Y[i_t] - 1) - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2 ) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
		        					    }
		        				      else { // obsT[k] < time[i_t]
		        					      c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to 
		        					    }
		            				}
											}
										else 		
											{//Reduced sample KM
												if (obsT[k] > time[i_t]){
		              				c_kl_int = D[i_t] / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) );
		            				}
		            				else {
		              				if (obsT[k] == time[i_t]){
		                				c_kl_int = (D[i_t] - status[k]) / ( ( Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2 ) * (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) );
		              				}
		              				else { // obsT[k] < time[i_t]
		                				c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
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
														na_k = D[i_t] / (Y[i_t] - 1 - C[i_t]/2);
														c_k_int = ( (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - 1 - D[i_t] - C[i_t]/2 ) );
													}
												else{
														if (obsT[l]==time[i_t])
															{
																// decrease the number of events if l was an event, and decrease the number at risk because l was in the risk set at time[i_t] 
																na_k = (D[i_t] - status[l]) / (Y[i_t] - 1 - (C[i_t] - (1-status[l]))/2);
																c_k_int = ( (Y[i_t] - 1 - (C[i_t] - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - 1 - ( D[i_t] - status[l] ) - (C[i_t] - (1-status[l]))/2 ) );  
															}
														else{ // obsT[l]<time[i_t]
																// do nothing since l is no longer in the risk set 
																na_k = D[i_t] / (Y[i_t] - C[i_t]/2);
																c_k_int = D[i_t] / ( (Y[i_t] - C[i_t]/2) * ( Y[i_t] - D[i_t] - C[i_t]/2 ) ); // ( (Y[i_t] - C[i_t]/2) * D[i_t]) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - D[i_t] - C[i_t]/2 ) ) reduces to
															}
													}
												Sk[l] = Sk[l] * (1 - na_k);
												covSSk[l] = covSSk[l] + c_k_int;
												}
											// k!= 0
											if (obsT[k] > time[i_t]){
								        if (obsT[l] > time[i_t]){
						          		c_kl_int = ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 2 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
						        		}
						        		else {
						          		if (obsT[l] == time[i_t]){
						          		  c_kl_int = ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 2 - (C[i_t] - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (D[i_t] - status[l]) - (C[i_t] - (1-status[l]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
						          		}
						          		else { // obsT[l] < time[i_t]
						          		  c_kl_int = ( (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to 
						          		}
						        		}
						      		}
						      		else {
						        		if (obsT[k] == time[i_t]){
						          		if (obsT[l] > time[i_t]){
						            		c_kl_int = ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 2 - (C[i_t] - (1-status[k]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2) * (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
						          		}
						          		else {
						            		if (obsT[l] == time[i_t]){
						              		c_kl_int = ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 2 - (C[i_t] - (1-status[k]) - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2) * (Y[i_t] - 1 - (D[i_t] - status[l]) - (C[i_t] - (1-status[l]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
						            		}
						            		else { // obsT[l] < time[i_t]
						              		c_kl_int = ( (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2) * (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to
						            		}
			 		                }	
								        }
								        else { // obsT[k] < time[i_t]
								          if (obsT[l] > time[i_t]){
								            c_kl_int = ( (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to 
								          }
								          else {
								            if (obsT[l] == time[i_t]){
								              c_kl_int = ( (Y[i_t] - 1 - (C[i_t] - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[l]) - (C[i_t] - (1-status[l]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (C[i_t] - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (D[i_t] - status[l]) - (C[i_t] - (1-status[l]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to
								            }
								            else { // obsT[l] < time[i_t]
								              c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to 
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
									S = S * (1 - D[i_t] / (Y[i_t] - C[i_t]/2));
									varS = varS + D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
									for (k=0;k<*N;k++){
										v1 = v1 + 1 / Sp[t*(*N) + k];
										if (k==0) {
											if (obsT[k]>time[i_t])
												{
													// decrease the number at risk because individual k was at risk at time[i_t] 
													na_k = D[i_t] / (Y[i_t] - 1 - C[i_t]/2);
													c_k_int = ( (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - 1 - D[i_t] - C[i_t]/2 ) );
												}
											else{
													if (obsT[k]==time[i_t])
														{
															// decrease the number of events if k was an event, and decrease the number at risk because k was in the risk set at time[i_t] 
															na_k = (D[i_t] - status[k]) / (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2);
															c_k_int = ( (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) * D[i_t] ) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - 1 - ( D[i_t] - status[k] ) - (C[i_t] - (1-status[k]))/2 ) );  
														}
													else{ // obsT[k]<time[i_t]
															// do nothing since k is no longer in the risk set 
															na_k = D[i_t] / (Y[i_t] - C[i_t]/2);
															c_k_int = D[i_t] / ( (Y[i_t] - C[i_t]/2) * ( Y[i_t] - D[i_t] - C[i_t]/2 ) ); // ( (Y[i_t] - C[i_t]/2) * D[i_t]) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - D[i_t] - C[i_t]/2 ) ) reduces to
														}
												}
											Sk[k] = Sk[k] * (1 - na_k);
											covSSk[k] = covSSk[k] + c_k_int;
											} // end of k==0 case
										// l == k case
										if (*varSkMethod == 1)
											{//precise
												if (obsT[k] > time[i_t]){
					 	              c_kl_int = (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) * D[i_t] / ( ( (Y[i_t] - 1) - D[i_t] - C[i_t]/2 ) * ( (Y[i_t] - 1) - D[i_t] - C[i_t]/2 ) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
		      					    }
		        				    else {
		        				      if (obsT[k] == time[i_t]){
		        					      c_kl_int = (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) * D[i_t] / ( ( (Y[i_t] - 1) - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2 ) * ( (Y[i_t] - 1) - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2 ) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
		        					    }
		        				      else { // obsT[k] < time[i_t]
		        					      c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to
		        					    }
		            				}
											}
										else 		
											{//Reduced sample KM
												if (obsT[k] > time[i_t]){
		              				c_kl_int = D[i_t] / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) );
		            				}
		            				else {
		              				if (obsT[k] == time[i_t]){
		                				c_kl_int = (D[i_t] - status[k]) / ( ( Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2 ) * (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) );
		              				}
		              				else { // obsT[k] < time[i_t]
		                				c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
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
														na_k = D[i_t] / (Y[i_t] - 1 - C[i_t]/2);
														c_k_int = ( (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - 1 - D[i_t] - C[i_t]/2 ) );
													}
												else{
														if (obsT[l]==time[i_t])
															{
																// decrease the number of events if l was an event, and decrease the number at risk because l was in the risk set at time[i_t] 
																na_k = (D[i_t] - status[l]) / (Y[i_t] - 1 - (C[i_t] - (1-status[l]))/2);
																c_k_int = ( (Y[i_t] - 1 - (C[i_t] - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - 1 - ( D[i_t] - status[l] ) - (C[i_t] - (1-status[l]))/2 ) );  
															}
														else{ // obsT[l]<time[i_t]
																// do nothing since l is no longer in the risk set 
																na_k = D[i_t] / (Y[i_t] - C[i_t]/2);
																c_k_int = D[i_t] / ( (Y[i_t] - C[i_t]/2) * ( Y[i_t] - D[i_t] - C[i_t]/2 ) ); // ( (Y[i_t] - C[i_t]/2) * D[i_t]) / ( (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * ( Y[i_t] - D[i_t] - C[i_t]/2 ) ) reduces to
															}
													}
												Sk[l] = Sk[l] * (1 - na_k);
												covSSk[l] = covSSk[l] + c_k_int;
												}
											// k!= 0
											if (obsT[k] > time[i_t]){
								        if (obsT[l] > time[i_t]){
						          		c_kl_int = ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 2 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
						        		}
						        		else {
						          		if (obsT[l] == time[i_t]){
						          		  c_kl_int = ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 2 - (C[i_t] - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (D[i_t] - status[l]) - (C[i_t] - (1-status[l]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
						          		}
						          		else { // obsT[l] < time[i_t]
						          		  c_kl_int = ( (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to 
						          		}
						        		}
						      		}
						      		else {
						        		if (obsT[k] == time[i_t]){
						          		if (obsT[l] > time[i_t]){
						            		c_kl_int = ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 2 - (C[i_t] - (1-status[k]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2) * (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
						          		}
						          		else {
						            		if (obsT[l] == time[i_t]){
						              		c_kl_int = ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 2 - (C[i_t] - (1-status[k]) - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2) * (Y[i_t] - 1 - (D[i_t] - status[l]) - (C[i_t] - (1-status[l]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) );
						            		}
						            		else { // obsT[l] < time[i_t]
						              		c_kl_int = ( (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (C[i_t] - (1-status[k]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[k]) - (C[i_t] - (1-status[k]))/2) * (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to
						            		}
			 		                }	
								        }
								        else { // obsT[k] < time[i_t]
								          if (obsT[l] > time[i_t]){
								            c_kl_int = ( (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to 
								          }
								          else {
								            if (obsT[l] == time[i_t]){
								              c_kl_int = ( (Y[i_t] - 1 - (C[i_t] - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - 1 - (D[i_t] - status[l]) - (C[i_t] - (1-status[l]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (C[i_t] - (1-status[l]))/2) * D[i_t] ) / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - 1 - (D[i_t] - status[l]) - (C[i_t] - (1-status[l]))/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to
								            }
								            else { // obsT[l] < time[i_t]
								              c_kl_int = D[i_t] / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ); // ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * D[i_t] ) / ( (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - D[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) * (Y[i_t] - C[i_t]/2) ) reduces to
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
				//printf("\n Sk[k] = ");printf("%f ",Sk[k]);
				//printf("\n S_Oi[k] = ");printf("%f \n",S_Oi[k]);
				for (l=k;l<*N;l++){
					//printf(" l = ");printf("%d ",l);
					printf("%f ",covSkSl[(*N)*(*N+1)/2 - (*N-k)*(*N+1-k)/2 + (l-k)]);			
					}		
			//printf("\n S_Oi[k] ");printf("%f ",S_Oi[k]);
			//printf("\n Sp[k] ");printf("%f ",Sp[t*(*N) + k]);
			//printf("\n S_Oi[k] / Sp[k] ");printf("%f ",S_Oi[k] / Sp[t*(*N) + k]);
			}*/
			S_E[t] = S_E[t] / (*N);
			//printf("\n varS je ");printf("%f ", varS);
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

// gcc my_source.c -o my_app
// ./my_app

#include <stdio.h>
int main()
{	
	double Y[] = {10, 8, 5, 4, 2};
	double D[] = {1, 3, 1, 1, 1};
	double C[] = {1, 0, 0, 1, 0};
	double time[] = { 300, 600, 900, 1200, 1500};
	double obsT[] = {300, 300, 600, 600, 600, 900, 1200, 1200, 1500, 3652};
	double status[] = {1, 0, 1, 1, 1, 1, 1, 0, 1, 0};
	int N = 10;
	int NT = 5;
	int neval = 6;
	int ind[] = {0, 1, 2, 3, 4, 4};
	double evalT[] = {300, 600, 900, 1200, 1500, 3652};
	double Sp[] = {0.9687181, 0.9800106, 0.9963772, 0.9972814, 0.9698040, 0.9837513, 0.9950476, 0.9748046, 0.9899375, 0.9053581, 0.9361604, 0.9570018, 0.9929095, 0.9945589, 0.9385132, 0.9672520, 0.9902263, 0.9514637, 0.9800204, 0.8125505, 0.9038967, 0.9338691, 0.9897334, 0.9918604, 0.9053689, 0.9475236, 0.9858984, 0.9282995, 0.9683776, 0.7171606, 0.8733717, 0.9128970, 0.9861086, 0.9890807, 0.8721212, 0.9257198, 0.9814509, 0.9041112, 0.9561258, 0.6296520, 0.8413325, 0.8898962, 0.9829739, 0.9863962, 0.8361351, 0.9032317, 0.9761966, 0.8784248, 0.9430540, 0.5507172, 0.5662581, 0.7174663, 0.9505583, 0.9595669, 0.4885878, 0.7114414, 0.9272684, 0.6915124, 0.8465829, 0.1478734};

	double S_E[neval];
	double varS_E[neval];
	int varSkMethod = 1; //"precise"
	//int varSkMethod = 2; //"reduced sample (KM)"
	double covSkSl[N*(N+1)/2];

	double *Y_p = &Y[0];
	double *D_p = &D[0]; 
	double *C_p = &C[0]; 	
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

	pseudoWH_pre2_LT(Y_p, D_p, C_p, time_p, obsT_p, status_p, S_E_p, varS_E_p, Sp_p, N_p, NT_p, neval_p, ind_p, varSkMethod_p, covSkSl_p);
	
	printf("\nPo evaluaciji.\n");	
	printf("Vrednost S_E je:\n");
	for (i=0;i<neval;i++){
		printf("%d ",i); 
		printf("%f ",S_E[i]);
	}
	printf("\n");

	return 0;
}
