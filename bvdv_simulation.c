#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <R.h>

//---------------------------------------------------------------------------------------------------
// Structure Setup
typedef struct animal_node {
	int demo_group;		// 0: Heifer, 1: Calf
	int age;			// Age (day old) of cattle.
	int days_to_cull;	// Days left for culling due to empty/test positiveness.
	int is_pregnant;	// 0: No, 1: Yes.
	int days_to_calv;	// Gestation period of 281 days.
	int days_to_estrous;// Estrous return every U(18, 24) days.
	int inf_stat;		// 0: Susceptible, 1: TI, 2: Recovered, 3: iMmuned, 4: PI
	int days_left;		// Days left; for the depletion of MAb U(112, 190), or to get recovered U(10, 20).
	int gest_at_inf;	// Gestation day at infection
	struct animal_node *next_node;
} animal_node;

//---------------------------------------------------------------------------------------------------
// MAIN FUNCTION
void bvdv_simulation(double *beta, double *rho, double *mu, int *tau, int *n_hfr, int *init_age, int *n_sample, int *f_psm, int *s_psm, int *mate_period, 
					int *d_test01, int *d_test02, int *n_sim_pos01, int *n_sim_pos02) {

	int _Random_binomial(int n, double p);
	double _Random_normal(double mu, double sd);
	double _Random_uniform(double a, double b);
	animal_node *_Add_animal_to_list(animal_node *new_node, animal_node *animal_list);
	void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);
	void _Create_animal_list(int n_farms, int **farm_info, int **group_size, double *mu, animal_node *animal_list[]);
	void _Seed_infection(int n_farms, int **farm_info, int **group_size, double *rho, animal_node *animal_list[]);
	void _Add_new_calf(int farm_id, int **group_size, int inf_stat, int gest_at_inf, animal_node *animal_list[]);
	void _Update_demo(int day, int n_farms, int **farm_info, int **group_size, animal_node *animal_list[]);
	int _Abortion(int inf_stat, int days_to_calv);
	void _Update_infection(int n_farms, int **farm_info, int **group_size, double *beta, animal_node *animal_list[]);
	void _Sample_serostatus(int farm_id, int *pnt_disease[][15], int *n_sample, animal_node *animal_list[]);
	void _Update_cull(int n_farms, int **group_size, animal_node *animal_list[]);
	
	srand(1);
	int i, j, day, n_farms= 1;
	double Se= 0.991, Sp= 0.813;
	int n_ptr= 0;

	// Array to store farm management information	
	int **farm_info;
	farm_info= calloc(n_farms, sizeof(int *));
	// Array to record group size for each farm
	int **group_size; // Group size for each farm
	group_size= calloc(n_farms, sizeof(int *));
	// Array of pointers of individual's infectious status
	int *pnt_disease[n_farms][15];
	// Array of test results of sampled animals
	int **arr_test;
	arr_test= calloc(n_farms, sizeof(int *));
	//int **arr_test[n_farms][(*n_sample)];
	// Array of first test result
	int *f_test;
	f_test= calloc(n_farms, sizeof(int));
	// Array of second test result
	int *s_test;
	s_test= calloc(n_farms, sizeof(int));

	for (i= 0; i< n_farms; i++) {
		farm_info[i]= calloc(5, sizeof(int));

		farm_info[i][0]= (*n_hfr); 		// Heifer herd size
		farm_info[i][1]= (*f_psm);		// Day of the first planned start date of mating (PSM)
		farm_info[i][2]= (*s_psm); 		// Day of the second PSM
		farm_info[i][3]= (*mate_period);// Day of weaning
		farm_info[i][4]= (*init_age); 	// Prospected age of heifers when they were weaned

		group_size[i]= calloc(10, sizeof(int)); // Hfr: S, T, R, M, P + C: S, T, R, M, P
		arr_test[i]= calloc((*n_sample), sizeof(int));

		for (j= 0; j< 15; j++) {
			pnt_disease[i][j]= &n_ptr;
		}		

	}

	// Declare an animal list
	animal_node *animal_list[n_farms]; 


	//********** Disease simulation **********//
	// Initialise demographics
	_Create_animal_list(n_farms, farm_info, group_size, mu, animal_list);
	//printf("Day= %i, Animals= S: %3i, T: %3i, R: %3i, M: %3i, P: %3i.\n", day, group_size[0][0], group_size[0][2], group_size[0][4], group_size[0][6], group_size[0][8]);
	
	// Simulate transmission
	for(day= 1; day <= (*d_test02); day++) {
		_Update_demo(day, n_farms, farm_info, group_size, animal_list);
		//printf("Day= %i, Animals= S: %3i, T: %3i, R: %3i, M: %3i, P: %3i.\n", day, group_size[0][0], group_size[0][2], group_size[0][4], group_size[0][6], group_size[0][8]);

		if (day == (*tau)) {
		    // Introduce BVDV
		    _Seed_infection(n_farms, farm_info, group_size, rho, animal_list);
			//printf("Day= %i, Animals= S: %3i, T: %3i, R: %3i, M: %3i, P: %3i.\n", day, group_size[0][0], group_size[0][2], group_size[0][4], group_size[0][6], group_size[0][8]);
		}

		if (day >= (*tau)) {
			_Update_infection(n_farms, farm_info, group_size, beta, animal_list);
			//printf("Day= %i, Animals= S: %3i, T: %3i, R: %3i, M: %3i, P: %3i.\n", day, group_size[0][0], group_size[0][2], group_size[0][4], group_size[0][6], group_size[0][8]);

			if (day == (*d_test01)) {
				for (i= 0; i< n_farms; i++) {
					_Sample_serostatus(i, pnt_disease, n_sample, animal_list);

					for (j= 0; j< (*n_sample); j++) {
						if ((*pnt_disease[i][j]) == 2) {
							arr_test[i][j]= _Random_binomial(1, Se);
						} else {
							arr_test[i][j]= _Random_binomial(1, 1.0-Sp);
						}
						f_test[i] += arr_test[i][j];
					}
				}
			}

			if (day == (*d_test02)) {
				for (i= 0; i< n_farms; i++) {
					for (j= 0; j< (*n_sample); j++) {
						if (arr_test[i][j] == 0) {
							if (*pnt_disease[i][j] == 2) {
								s_test[i] += _Random_binomial(1, Se);
							} else {
								s_test[i] += _Random_binomial(1, 1.0-Sp);
							}
						}
					}
				}
			}

			_Update_cull(n_farms, group_size, animal_list);
		} else {
			_Update_cull(n_farms, group_size, animal_list);
		}
	}


	//********** Update the results **********//
	(*n_sim_pos01)= f_test[0];
	(*n_sim_pos02)= s_test[0];

	
	//********** Clear the memory **********//
	for (i= 0; i< n_farms; i++) {
		animal_node *eraser; // Declare a pointer to free the animal list
		eraser= animal_list[i];
		while(animal_list[i]) {
			eraser= animal_list[i];
			animal_list[i]= eraser->next_node;
			free(eraser);
		}
		free(farm_info[i]);
		free(group_size[i]);
	}
	free(farm_info);
	free(group_size);
	free(arr_test);
	free(f_test);
	free(s_test);

}

// Number of success from binomial distribution: It retunrs a number of success according to a probability "p".
int _Random_binomial(int n, double p) {
    int index;
	int n_pos= 0; 
	int i= 0;
	if (p > 1) {p= 1.0;}
	if (p < 0) {p= 0.0;}
	
	double x;

	if (p == 0.0) {
		return 0;
	} else if (p == 1.0) {
		return n;
	} else {
	    do {
	       x= (double)rand()/(double)RAND_MAX;
	       if(p >= x) {index= 1;} else {index= 0;}
	       n_pos += index;
	       i++;
	    } while (i < n);
	
		return n_pos;
	}
}

// Random number from normal distribution: It retunrs a random number according to a normal distribution ("mu", "sigma").
// CAUTION!!!!!!!!!!!!!!!!!!!! SOMETHING IS WRONG WITH THIS NORMAL DISTRIBUTION GENERATOR.
// PREV VERSION: I allowed "u" being zero, so the returned value can be -Inf.
// I fixed it, but need to test it before use it.
double _Random_normal(double mu, double sd) {
	double two_pi= 2.0 * acos(-1.0); // acos(-1) is PI
	double u, v, z;

	do {
		u= (double)rand()/(double)RAND_MAX;
		v= (double)rand()/(double)RAND_MAX;
	} while (u == 0.0);

	z= sqrt(-2.0 * log(u)) * cos(two_pi * v);

	return (mu + z * sd);
}

// Random number between two numbers: It retunrs a random number according to a uniform distribution ("a", "b").
double _Random_uniform(double a, double b) {
	double small;
	double diff= fabs(b - a);
	double x;

	if (a > b) {small= b;} else if (a < b) {small= a;} else {return a;}

	x= diff * (double)rand() / (double)(RAND_MAX/1.0);
	small= small + x;

	return small;
}

// Add new animal to the list by age: It append a node of animal to the existing list of animals according to the animal's age.
animal_node *_Add_animal_to_list(animal_node *new_node, animal_node *animal_list) {
	animal_node *previous_node;
	animal_node *current_node;

	previous_node= animal_list;
	current_node= animal_list;

	if (current_node == NULL) {
		animal_list= new_node;
		return (animal_list);
	} else if (current_node->demo_group == -99999) {
		new_node->next_node= current_node;
		animal_list= new_node;
		return (animal_list);
	} else {
		while (current_node->demo_group == 99999) {
			previous_node= current_node;
			current_node= current_node->next_node;
		}
		previous_node->next_node= new_node;
		new_node->next_node= current_node;
		return (animal_list);
	}
}

// Change group size based on the demographic group: It changes the number of animals in a certain group (as a combination of "demo_group" and "inf_stat").
void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals) {
	group_size[farm_id][(inf_stat * 2 + demo_group)] += n_animals;
}

// Create animal list: It generates a list of cattle for inital demographic setting of simulation.
void _Create_animal_list(int n_farms, int **farm_info, int **group_size, double *mu, animal_node *animal_list[]) {     
	int _Random_binomial(int n, double p);
	animal_node *_Add_animal_to_list(animal_node *new_node, animal_node *animal_list);
    void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);  

	int i, j;
	animal_node *new_node;

	for (i= 0; i< n_farms; i++) {
		new_node= (animal_node *)malloc(sizeof(animal_node));
		new_node->demo_group= -99999;
		new_node->next_node= NULL;
		animal_list[i]= new_node;
	
		new_node= (animal_node *)malloc(sizeof(animal_node));
	    new_node->demo_group= 99999;
	    new_node->next_node= NULL;
		animal_list[i]= _Add_animal_to_list(new_node, animal_list[i]);

		// Young heifers
        for (j= 0; j< farm_info[i][0]; j++) {
			new_node= (animal_node *)malloc(sizeof(animal_node));
	
			new_node->demo_group= 0;
			new_node->age= farm_info[i][4];
			new_node->days_to_cull= 99999;
			new_node->is_pregnant= 0;
			new_node->days_to_calv= 0;
			new_node->days_to_estrous= 0;
			if (_Random_binomial(1, (*mu)) == 1) {
				new_node->inf_stat= 2;
			} else {
				new_node->inf_stat= 0;
			}
			new_node->days_left= 0;
			new_node->gest_at_inf= 0;
		    new_node->next_node= NULL;
			animal_list[i]= _Add_animal_to_list(new_node, animal_list[i]);
			_Change_group_size(group_size, i, new_node->demo_group, new_node->inf_stat, 1);  
		}
	}
}

// Infection seed: It initiates BVD by randomly converting X number of animals into PI.
void _Seed_infection(int n_farms, int **farm_info, int **group_size, double *rho, animal_node *animal_list[]) {
    void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);  

	int i, n_pi;

	animal_node *previous_node;
	animal_node *current_node;

	for (i= 0; i< n_farms; i++) {
		previous_node= animal_list[i];
		current_node= animal_list[i];

		n_pi= (int)round((double)farm_info[i][0] * (*rho));

		while (current_node != NULL) {
			if (current_node->demo_group == 99999 || current_node->demo_group == -99999) {
				previous_node= current_node;
				current_node= current_node->next_node;
			} else if (current_node->demo_group == 0 && current_node->inf_stat == 0) {
				if (n_pi > 0) {
					_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
					current_node->inf_stat= 4;
					n_pi -=1;
					_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
					previous_node= current_node;
					current_node= current_node->next_node;
				} else {
					previous_node= current_node;
					current_node= current_node->next_node;
				}
			} else {
				previous_node= current_node;
				current_node= current_node->next_node;
			}
		}
	}
}

// Create a new calf: It appends a node of new calf to a list of cattle according to its dam's status.
void _Add_new_calf(int farm_id, int **group_size, int inf_stat, int gest_at_inf, animal_node *animal_list[]) {
	int _Random_binomial(int n, double p);
	double _Random_normal(double mu, double sd);
	animal_node *_Add_animal_to_list(animal_node *new_node, animal_node *animal_list);
    void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);  

	animal_node *new_node;
	
	new_node= (animal_node *)malloc(sizeof(animal_node));

	new_node->demo_group= 1;
	new_node->age= 0;
	new_node->days_to_cull= 99999;
	new_node->is_pregnant= 0;
	new_node->days_to_calv= 0;
	new_node->days_to_estrous= 0;
	new_node->gest_at_inf= 0;
    new_node->next_node= NULL;

	// Define infection status
	// 0: Susceptible, 1: TI, 2: Recovered, 3: iMmuned, 4: PI
	// From Susceptible: Susceptible calves
	if (inf_stat == 0) {
		new_node->inf_stat= 0;
	// From PI: PI calves (abortion rate was assumed to be the same across TI and PI)
	} else if (inf_stat == 4) {
		new_node->inf_stat= 4;
	// From Infected, Recovered or immuned: Infection status of calves depends on the time of infection
	} else {
		// From recovered or vaccinated (Infection status cannot be used as a dam can be infected during risk period but recovered by the time of delivery).
		// If dams were recovered or vaccinated but never infected, gest_at_inf should be 0.
		if (gest_at_inf == 0) {
			new_node->inf_stat= 3; // When dams are recovered, new calves will get maternal Ab
			new_node->days_left= (int)round(_Random_normal(155.0, 31.0));
		} else {
			// Early pregnancy
			if (gest_at_inf >= 239) { // Day 0 ~ Day 41
				new_node->inf_stat= 3; // Immuned calf from TI cattle
				new_node->days_left= (int)round(_Random_normal(155.0, 31.0));
			// Mid pregnancy
			} else if (gest_at_inf >= 130 && gest_at_inf < 239) { // Day 42 ~ Day 150 (Pregnancy period was categorised as same as Damman et al., 2015)
				if (_Random_binomial(1, 0.933) == 1) { // 0.933 is conditional probability of calving PI (Cond: No abortion & CD)
					new_node->inf_stat= 4;
				} else {
					if (_Random_binomial(1, 0.5) == 1) { // 0.50 is conditional probability of calving immuned (Cond: No abortion, No PI)
						new_node->inf_stat= 3; // Immuned calf from TI cattle
						new_node->days_left= (int)round(_Random_normal(155.0, 31.0));
					} else {
						new_node->inf_stat= 2; // Recovered calf from TI cattle
					}
				}
			// Late pregnancy
			} else {
				new_node->inf_stat= 2;
			}
		}
	}

	_Change_group_size(group_size, farm_id, new_node->demo_group, new_node->inf_stat, 1);
	animal_list[farm_id]= _Add_animal_to_list(new_node, animal_list[farm_id]);
}

// Update demographic status: It updates animal's demographic compartment for every farm according to farm's demographic event date. 
void _Update_demo(int day, int n_farms, int **farm_info, int **group_size, animal_node *animal_list[]) {
    void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);  
	int _Random_binomial(int n, double p);
	double _Random_uniform(double a, double b);
	void _Add_new_calf(int farm_id, int **group_size, int inf_stat, int gest_at_inf, animal_node *animal_list[]);

	int i;
	int f_psm, s_psm, mate_period;
	double p_calf_mort_hfr, p_conc_bull, rr_conc;

	p_calf_mort_hfr= 0.00086;
	rr_conc= 0.6805;

	animal_node *previous_node;
	animal_node *current_node;

	for (i= 0; i< n_farms; i++) {
		f_psm= farm_info[i][1];
		s_psm= farm_info[i][2];
		mate_period= farm_info[i][3];

		previous_node= animal_list[i];
		current_node= animal_list[i];

		p_conc_bull= 0.61;
		
		while (current_node != NULL) {
			if (current_node->demo_group == 99999 || current_node->demo_group == -99999) {
				previous_node= current_node;
				current_node= current_node->next_node;
			} else {
				// Calves
				if (current_node->demo_group == 1) {
					current_node->days_to_cull -=1;
					current_node->age +=1;

					// Natural mortality of until weaning
					if (_Random_binomial(1, p_calf_mort_hfr) == 1) {
						current_node->days_to_cull= 0;
					}

					// Move to the next node
					previous_node= current_node;
					current_node= current_node->next_node;
				// Young heifers
				} else if (current_node->demo_group == 0) {
					current_node->days_to_cull -=1;
					current_node->age +=1;

					// At the puberty
					if (current_node->age == 380) { // Puberty starts at 380 days old. McNaughton et al., 2002
						current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // Set the counter upto the second heat (silent first heat)
					// After the puberty
					} else if (current_node->age > 380) {
						// Pregnant heifers
						if (current_node->is_pregnant == 1) {
							// Natural abortion
							if (_Random_binomial(1, 0.00013) == 1) {
								current_node->is_pregnant= 0;
								current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0));// Set the counter upto the second heat (womb preparation + silent first heat)
								current_node->days_to_calv= 0;
								current_node->gest_at_inf= 0;
							} else {
								// Calving
								if (current_node->days_to_calv == 0) {
									_Add_new_calf(i, group_size, current_node->inf_stat, current_node->gest_at_inf, animal_list);
									current_node->is_pregnant= 0;
									current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0)); 
									// Set the counter upto the second heat (womb preparation (calving to first oestrous (Olds et al., 1953) + silent first heat)
									current_node->gest_at_inf= 0;
								} else {
									current_node->days_to_calv -=1; // REDUCE THE DAYS TO CALV
								}
							}
						// Empty heifers
						} else {
							// On Heat
							if (current_node->days_to_estrous == 0) {
								// DURING THE 9 WEEKS OF BREEDING PERIOD
								if ((day >= f_psm && day < f_psm + mate_period) || (day >= s_psm && day < s_psm + mate_period)) {
									// TI or PI heifers
									if (current_node->inf_stat == 1 || current_node->inf_stat == 4) {
										if (_Random_binomial(1, p_conc_bull * rr_conc) == 1) {
											current_node->is_pregnant= 1;
											current_node->days_to_calv= 280; // It means the calving will occurr after 281 days as calving occurs under the condition of days_to_calv = 0.
										} else {
											current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // RESET THE DAYS TO ESTROUS DEPENDING ON THE PREGNANT STATUS
										}
									} else {
										if (_Random_binomial(1, p_conc_bull) == 1) {
											//printf("Conc.\n");
											current_node->is_pregnant= 1;
											current_node->days_to_calv= 280; // It means the calving will occurr after 281 days as calving occurs under the condition of days_to_calv = 0.
										} else {
											current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // RESET THE DAYS TO ESTROUS DEPENDING ON THE PREGNANT STATUS
										}
									}
								// Out of mating period
								} else {
									current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // RESET THE DAYS TO ESTROUS DEPENDING ON THE PREGNANT STATUS
								}
							// No heat
							} else {
								current_node->days_to_estrous -=1; // REDUCE THE DAYS TO ESTROUS
							}
						}
					}
					
					// Move to the next node
					previous_node= current_node;
					current_node= current_node->next_node;
				}
			}
		}
	}
}

// Abortion: It returns an abortion indicater according to the gestation period (= "days_to_calv") of infected cattle. We assumed that there was no difference of abortion rate between TI and PI.
int _Abortion(int inf_stat, int days_to_calv) {
	int _Random_binomial(int n, double p);

	// TI
	if (inf_stat == 1) {
		// Early pregnancy (abortion probability= 0.1212)
		if (days_to_calv >= 239) {
			if (_Random_binomial(1, 0.1212) == 1) {return 1;
			} else {return 0;}
		// Mid pregnancy (abortion probability= 0.1194, congenital deformality= 0.0625)
		} else if (days_to_calv >= 130 && days_to_calv < 239) {
			if (_Random_binomial(1, 0.1744) == 1) {return 1;
			} else {return 0;}
		// Late pregnancy: No abortion
		} else {return 0;}
	// PI
	} else if (inf_stat == 4) {
		// Early pregnancy (abortion probability= 0.1212: abortion probability/day= 0.0086)
		if (days_to_calv >= 239) {
			if (_Random_binomial(1, 0.0086) == 1) {return 1;
			} else {return 0;}
		// Mid pregnancy (abortion probability= 0.1194, congenital deformality= 0.0625: abortion+CD probability/day= 0.0128)
		} else if (days_to_calv >= 130 && days_to_calv < 239) {
			if (_Random_binomial(1, 0.0128) == 1) {return 1;
			} else {return 0;}
		// Late pregnancy: No abortion
		} else {return 0;}
	} else {return 0;
	}
}

// Disease spread: It updates BVD infection status of each cattle.
void _Update_infection(int n_farms, int **farm_info, int **group_size, double *beta, animal_node *animal_list[]) {
	// 0: Calves (F and M), 1: RHs1 (weaning to off-grazing), 2: RHs2 (separately managed), 3: MAs, 4: R1 bulls, 5: R2 bulls, 6: R3 bulls.
	int _Random_binomial(int n, double p);
	double _Random_uniform(double a, double b);
    void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);  
	int _Abortion(int inf_stat, int days_to_calv);
	
	// Variables for simulation
	double beta_ti= (*beta) * 0.05; // Low= 1%, Mod= 5%, High= 10%

	int i, j, n_herd, n_pi, n_ti;
	double prev_pi, prev_ti, p_inf;
	double lambda_pi= 0.0019;
	
	animal_node *previous_node;
	animal_node *current_node;

	for (i= 0; i< n_farms; i++) {
		n_herd= 0;
		n_pi= 0;
		n_ti= 0;

		for (j= 0; j< 10; j++) {
			n_herd += group_size[i][j];
		}

		n_pi= group_size[i][8] + group_size[i][9];
		n_ti= group_size[i][2] + group_size[i][3];
			
		prev_pi= (double)n_pi / (double)n_herd;
		prev_ti= (double)n_ti / (double)n_herd;
			
		previous_node= animal_list[i];
		current_node= animal_list[i];
	
		// Probability of infection per demographic group
		p_inf= 1 - exp(-1.0 * ((*beta)*prev_pi + beta_ti*prev_ti));
		
		// Update transmission
		while (current_node != NULL) {
			// 0: Susceptible, 1: TI, 2: Recovered, 3: iMmuned, 4: PI.
			// Recovered cattle
			if (current_node->demo_group == 99999 || current_node->demo_group == -99999 || current_node->inf_stat == 2) {
				previous_node= current_node;
				current_node= current_node->next_node;
			} else {
				// Susceptible -> TI
				if (current_node->inf_stat == 0) {
					if (_Random_binomial(1, p_inf) == 1) {
						// Update infection status
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						current_node->inf_stat= 1;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
						current_node->days_left= (int)round(_Random_uniform(10.0, 20.0));
						// Reproduction outcome
						if (current_node->is_pregnant == 1) {
							if (_Abortion(current_node->inf_stat, current_node->days_to_calv) == 1) {
								current_node->is_pregnant= 0;
								current_node->days_to_calv= 0;
								current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0)); 
								// Set the counter upto the second heat (womb preparation (calving to first oestrous (Olds et al., 1953) + silent first heat)
								current_node->gest_at_inf= 0;
							} else {
								current_node->gest_at_inf= current_node->days_to_calv;
							}
						}
					}
					previous_node= current_node;
					current_node= current_node->next_node;
				// TI -> Recovered
				} else if (current_node->inf_stat == 1) {
					// Recovery from TI
					if (current_node->days_left == 0) {
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						current_node->inf_stat= 2;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
					// Not yet recovered
					} else {
						current_node->days_left -=1;
					}
					previous_node= current_node;
					current_node= current_node->next_node;
				// Waning-off of Maternal Ab
				} else if (current_node->inf_stat == 3) {
					// iMmuned -> Susceptible
					if (current_node->days_left == 0) {
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						current_node->inf_stat= 0;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
					} else if (current_node->days_left > 0 && current_node->age >= 14) {
						// Anamnestic response
						if (_Random_binomial(1, p_inf) == 1) {
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
							current_node->inf_stat= 2;
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
							current_node->days_left= 0;
						} else {
							// Depletion of MAb
							current_node->days_left -=1;
						}
					} else {
						// Depletion of MAb
						current_node->days_left -=1;
					}
					previous_node= current_node;
					current_node= current_node->next_node;
				// Mortality of PIs
				} else if (current_node->inf_stat == 4) {
					// Death of PI
					if (_Random_binomial(1, 1 - lambda_pi) == 0) {
						current_node->days_to_cull= 0;
					// PI survived
					} else {
						// Reproduction outcome
						if (current_node->is_pregnant == 1 && _Abortion(current_node->inf_stat, current_node->days_to_calv) == 1) {
							current_node->is_pregnant= 0;
							current_node->days_to_calv= 0;
							current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0));
							// Set the counter upto the second heat (womb preparation (calving to first oestrous (Olds et al., 1953) + silent first heat)
							current_node->gest_at_inf= 0;
						}
					}
					previous_node= current_node;
					current_node= current_node->next_node;
				}
			}
		}
	}
}

// Cull/sell cattle based on the culling status
void _Update_cull(int n_farms, int **group_size, animal_node *animal_list[]) {
	int i;

	animal_node *previous_node;
	animal_node *current_node;

	for (i= 0; i< n_farms; i++) {
		previous_node= animal_list[i];
		current_node= animal_list[i];

		while (current_node != NULL) {
			if (current_node->demo_group == 99999 || current_node->demo_group == -99999) {
				previous_node= current_node;
				current_node= current_node->next_node;
			} else {
				if (current_node->days_to_cull == 0) {
					previous_node->next_node= current_node->next_node;
					_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
					free(current_node);
					current_node= previous_node->next_node;
				} else {
					previous_node= current_node;
					current_node= current_node->next_node;
				}
			}
		}
	}
}

// Sample the heifers for the test and run ELISA
void _Sample_serostatus(int farm_id, int *pnt_disease[][15], int *n_sample, animal_node *animal_list[]) {
	int i= 0;

	animal_node *previous_node;
	animal_node *current_node;
	previous_node= animal_list[farm_id];
	current_node= animal_list[farm_id];

	while (current_node != NULL) {
		if (current_node->demo_group == 99999 || current_node->demo_group == -99999) {
			previous_node= current_node;
			current_node= current_node->next_node;
		} else {
//			if (current_node->demo_group == 0 && current_node->inf_stat != 4) {
			if (current_node->demo_group == 0) {
				if (i < (*n_sample)) {
					pnt_disease[farm_id][i]= &(current_node->inf_stat);
					i +=1;
				}
				previous_node= current_node;
				current_node= current_node->next_node;
			} else {
				previous_node= current_node;
				current_node= current_node->next_node;
			}
		}
	}
}
