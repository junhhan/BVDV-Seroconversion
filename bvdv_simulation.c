#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <R.h>

//---------------------------------------------------------------------------------------------------
// Structure Setup
typedef struct animal_node {
	int demo_group;		// 0: Heifer, 1: Calf
	int preg;			// 0: Negative, 1: Positive
	int days_to_calv;	// Days left for calving (calving at 0). Negative value is possible and it indicates that the cattle is empty.
	int inf_stat;		// 0: Susceptible, 1: TI, 2: Recovered, 3: iMmuned, 4: PI
	int days_to_sr;		// Days left to get suscep/recovered
	int gest_at_doi;	// Gestation day at the day of infection
	struct animal_node *next_node;
} animal_node;

//---------------------------------------------------------------------------------------------------
// MAIN FUNCTION
void bvdv_simulation(double *beta, double *rho, int *tau, int *n_hfr, int *n_sample, double *p_conc, int *d_mate_start, int *d_mate_end, 
					int *d_test01, int *d_test02, int *n_sim_pos01, int *n_sim_pos02, double *p_sim01, double *p_sim02) {
	double _Random_normal(double mu, double sd);
	int _Random_binomial(int n, double p);
	double _Random_uniform(double a, double b);
	animal_node *_Add_animal_to_list(animal_node *new_node, animal_node *animal_list);
	void _Add_group_size(int **group_size, int farm_id, int demo_group, int inf_stat);
	void _Draw_group_size(int **group_size, int farm_id, int demo_group, int inf_stat);
	void _Add_init_animal(int **group_size, int farm_id, animal_node *animal_list[]);
	void _Create_animal_list(int n_farms, int **group_size, animal_node *animal_list[], int *n_hfr);
	void _Seed_infection(int n_farms, double *rho, int **group_size, animal_node *animal_list[], int *n_hfr);
	void _Add_new_calf(int **group_size, int farm_id, animal_node *animal_list[], int inf_stat, int gest_at_doi);
	void _Update_demo_stat(int day, int n_farms, int **group_size, animal_node *animal_list[], int *d_mate_start, int *d_mate_end, double *p_conc);
	int _Abortion(int inf_stat, int days_to_calv);
	void _Update_inf_stat(double *beta, int n_farms, int **group_size, animal_node *animal_list[]);
	void _Sample_serostatus(animal_node *animal_list[], int *pnt_disease[], int *n_sample);
	

	int i, n_farms= 1;
	srand(time(NULL));

	int **group_size; // Group size for each farm
	group_size= calloc(n_farms, sizeof(int *));
	for (i= 0; i< n_farms; i++) {
		group_size[i]= calloc(10, sizeof(int)); // Hfr: S, T, R, M, P + C: S, T, R, M, P
	}

	int day;
	int *pnt_disease[(*n_sample)];
	int arr_test[(*n_sample)];
	double Se= 0.969, Sp= 0.978;

	animal_node *animal_list[n_farms]; // Declare an animal list


	//********** Disease simulation **********//
	// Initialise demographics
	_Create_animal_list(n_farms, group_size, animal_list, n_hfr); // Set initial values for animal lists
	//printf("Day= %i, Animals= S: %3i, T: %3i, R: %3i, M: %3i, P: %3i.\n", day, group_size[0], group_size[1], group_size[2], group_size[3], group_size[4]);
	
	// Simulate transmission
	for(day= 0; day< (*d_test02); day++) {
		if (day == (*tau)) {
		    // Introduce BVDV
			_Seed_infection(n_farms, rho, group_size, animal_list, n_hfr);
			//printf("Day= %i, Animals= S: %3i, T: %3i, R: %3i, M: %3i, P: %3i.\n", day, group_size[0], group_size[1], group_size[2], group_size[3], group_size[4]);
		}

		_Update_demo_stat(day, n_farms, group_size, animal_list, d_mate_start, d_mate_end, p_conc);
		//printf("Day= %i, Animals= S: %3i, T: %3i, R: %3i, M: %3i, P: %3i.\n", day, group_size[0], group_size[1], group_size[2], group_size[3], group_size[4]);

		if (day >= (*tau)) {
			_Update_inf_stat(beta, n_farms, group_size, animal_list);
			//printf("Day= %i, Animals= S: %3i, T: %3i, R: %3i, M: %3i, P: %3i.\n", day, group_size[0], group_size[1], group_size[2], group_size[3], group_size[4]);

			if (day == (*d_test01) - 1) {
				_Sample_serostatus(animal_list, pnt_disease, n_sample);

				for (i= 0; i< (*n_sample); i++) {
					//if (*pnt_disease[i] == 2) {arr_test[i]= 1;} else {arr_test[i]= 0;}
					//(*n_sim_pos01) += arr_test[i];
					(*p_sim01)= (double)group_size[0][2] / (double)(group_size[0][0] + group_size[0][1] + group_size[0][2] + group_size[0][3] + group_size[0][4] + 
								group_size[0][5] + group_size[0][6] + group_size[0][7] + group_size[0][8] + group_size[0][9]);

					if (*pnt_disease[i] == 2) {
						arr_test[i]= _Random_binomial(1, Se);
					} else {
						arr_test[i]= _Random_binomial(1, 1.0-Sp);
					}
					(*n_sim_pos01) += arr_test[i];
				}
			}

			if (day == (*d_test02) - 1) {
				(*n_sim_pos02)= (*n_sim_pos01);

				for (i= 0; i< (*n_sample); i++) {
					if (arr_test[i] == 0) {
						//if (*pnt_disease[i] == 2) {(*n_sim_pos02) +=1;}

						if (*pnt_disease[i] == 2) {
							(*n_sim_pos02) += _Random_binomial(1, Se);
						} else {
							(*n_sim_pos02) += _Random_binomial(1, 1.0-Sp);
						}
					} 
				}
				(*p_sim02)= (double)group_size[0][2] / (double)(group_size[0][0] + group_size[0][1] + group_size[0][2] + group_size[0][3] + group_size[0][4] + 
							group_size[0][5] + group_size[0][6] + group_size[0][7] + group_size[0][8] + group_size[0][9]);
			}
		}
	}
	
	for (i= 0; i< n_farms; i++) {
		animal_node *eraser; // Declare a pointer to free the animal list
		eraser= animal_list[i];
		while(animal_list[i]) {
			eraser= animal_list[i];
			animal_list[i]= eraser->next_node;
			free(eraser);
		}
	}

}

// Random from normal distribution (modified Box-Muller method)
double _Random_normal(double mu, double sd) {
	double two_pi= 2.0 * acos(-1); // acos(-1) is PI
	double u, v, z;

	u= (double)rand()/(double)RAND_MAX;
	v= (double)rand()/(double)RAND_MAX;

	z= sqrt(-2.0 * log(u)) * cos(two_pi * v);

	return (mu + z * sd);
}

// Random from binomial distribution
int _Random_binomial(int n, double p) {
    int index;
	int n_pos= 0; 
	int i= 0;
	
	double x;

    do {
       x= (double)rand()/(double)RAND_MAX;
       if(p >= x) {index= 1;} else {index= 0;}
       n_pos += index;
       i++;
    } while (i < n);

	return n_pos;
}

// Random from uniform distribution
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

// Add group size based on the demographic group: It increases the number of animals in a certain group (as a combination of "demo_group" and "inf_stat").
void _Add_group_size(int **group_size, int farm_id, int demo_group, int inf_stat) {
	group_size[farm_id][(demo_group * 5) + inf_stat] +=1;
}

// Subtract group size based on the demographic group: It decreases the number of animals in a certain group (as a combination of "demo_group" and "inf_stat").
void _Draw_group_size(int **group_size, int farm_id, int demo_group, int inf_stat) {
	group_size[farm_id][(demo_group * 5) + inf_stat] -=1;
}

// Generate an animal
void _Add_init_animal(int **group_size, int farm_id, animal_node *animal_list[]) {
	animal_node *_Add_animal_to_list(animal_node *new_node, animal_node *animal_list);
	void _Add_group_size(int **group_size, int farm_id, int demo_group, int inf_stat);
	
	animal_node *new_animal;
	new_animal= (animal_node *)malloc(sizeof(animal_node));

	new_animal->demo_group= 0;
	new_animal->preg= 0;
	new_animal->days_to_calv= 0;
	new_animal->inf_stat= 0;
	new_animal->days_to_sr= 0;
	new_animal->gest_at_doi= 0;
	new_animal->next_node= NULL;

	animal_list[farm_id]= _Add_animal_to_list(new_animal, animal_list[farm_id]);
	_Add_group_size(group_size, farm_id, new_animal->demo_group, new_animal->inf_stat);
}

// Create an animal list: It generates a list of cattle for inital demographic setting of simulation.
void _Create_animal_list(int n_farms, int **group_size, animal_node *animal_list[], int *n_hfr) {
	void _Add_init_animal(int **group_size, int farm_id, animal_node *animal_list[]);

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
	
	    for (j= 0; j< (*n_hfr); j++) {
	        _Add_init_animal(group_size, i, animal_list);
		}
	}
}

// Infection seed: It initiates BVD by randomly designating 10% of new born calves as PI.
void _Seed_infection(int n_farms, double *rho, int **group_size, animal_node *animal_list[], int *n_hfr) {
	animal_node *_Add_animal_to_list(animal_node *new_node, animal_node *animal_list);
	void _Add_group_size(int **group_size, int farm_id, int demo_group, int inf_stat);

	int i, j, n_pi;

	animal_node *new_animal;
	n_pi= (int)round((double)(*n_hfr) * (*rho));

	for (i= 0; i< n_farms; i++) {

	    for (j= 0; j< n_pi; j++) {
			new_animal= (animal_node *)malloc(sizeof(animal_node));
		
			new_animal->demo_group= 0;
			new_animal->preg= 0;
			new_animal->days_to_calv= 0;
			new_animal->inf_stat= 4;
			new_animal->days_to_sr= 0;
			new_animal->gest_at_doi= 0;
			new_animal->next_node= NULL;
	
			animal_list[i]= _Add_animal_to_list(new_animal, animal_list[i]);
			_Add_group_size(group_size, i, new_animal->demo_group, new_animal->inf_stat);
		}
	}	
}

// Create a new calf: It appends a node of new calf to a list of cattle according to its dam's status.
void _Add_new_calf(int **group_size, int farm_id, animal_node *animal_list[], int inf_stat, int gest_at_doi) {
	double _Random_uniform(double a, double b);
	int _Random_binomial(int n, double p);
	animal_node *_Add_animal_to_list(animal_node *new_node, animal_node *animal_list);
	void _Add_group_size(int **group_size, int farm_id, int demo_group, int inf_stat);

	animal_node *new_calf;
	new_calf= (animal_node *)malloc(sizeof(animal_node));

	new_calf->demo_group= 1;
	new_calf->preg= 0;
	new_calf->days_to_calv= 0;
	new_calf->gest_at_doi= 0;
	new_calf->next_node= NULL;

	// Define infection status
	// 0: Susceptible, 1: TI, 2: Recovered, 3: iMmuned, 4: PI
	// From Susceptible: Susceptible calves
	if (inf_stat == 0) {
		new_calf->inf_stat= 0;
	// From PI: PI calves (abortion rate was assumed to be the same across TI and PI)
	} else if (inf_stat == 4) {
		new_calf->inf_stat= 4;
	// From recovered
	} else if (inf_stat == 2 && gest_at_doi == 0) { // Dams had been Infected then Recovered by the time of calving, so adjust it by "gest_at_doi == 0")
		new_calf->inf_stat= 3;
		new_calf->days_to_sr= (int)round(_Random_uniform(150.0, 180.0));
	// From Infected: Infection status of calves depends on the time of infection
	} else {
		// Early pregnancy
		if (gest_at_doi >= 238) {
			new_calf->inf_stat= 3; // Immuned calf from TI cattle
			new_calf->days_to_sr= (int)round(_Random_uniform(150.0, 180.0));
		// Mid pregnancy
		} else if (gest_at_doi >= 129 && gest_at_doi < 238) {
			if ((int)_Random_binomial(1.0, 0.934) == 1) { // 0.934 is conditional probability of calving PI (Cond: No abortion)
				new_calf->inf_stat= 4;
			} else {
				if ((int)_Random_binomial(1.0, 0.5) == 1) { // 0.5 is conditional probability of calving immuned (Cond: No abortion, No PI)
					new_calf->inf_stat= 3; // Immuned calf from TI cattle
					new_calf->days_to_sr= (int)round(_Random_uniform(150.0, 180.0));
				} else {
					new_calf->inf_stat= 2; // Recovered calf from TI cattle
				}
			}
		// Late pregnancy
		} else {
			new_calf->inf_stat= 2;
		}
	}
	
	animal_list[farm_id]= _Add_animal_to_list(new_calf, animal_list[farm_id]);
	_Add_group_size(group_size, farm_id, new_calf->demo_group, new_calf->inf_stat);
}

// Update demographic status: It updates animal's demographic compartment for every farm according to farm's demographic event date.
void _Update_demo_stat(int day, int n_farms, int **group_size, animal_node *animal_list[], int *d_mate_start, int *d_mate_end, double *p_conc) {
	int _Random_binomial(int n, double p);
	void _Add_new_calf(int **group_size, int farm_id, animal_node *animal_list[], int inf_stat, int gest_at_doi);
	void _Add_group_size(int **group_size, int farm_id, int demo_group, int inf_stat);

	// Variables for simulation
	int i;
	double p_wean= 0.9;

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
				if (current_node->demo_group == 0) {
					if (current_node->preg == 1) {current_node->days_to_calv -=1;}
					// MATING
					if ((day >= *d_mate_start && day <= *d_mate_end) && current_node->preg == 0) {
						if (_Random_binomial(1, (*p_conc)) == 1) {
							current_node->days_to_calv= 279;
							current_node->preg= 1;
						}
					// CALVING
					} else if (current_node->days_to_calv == 0 && current_node->preg == 1) {
						if (_Random_binomial(1, p_wean) == 1) {_Add_new_calf(group_size, i, animal_list, current_node->inf_stat, current_node->gest_at_doi);}
						current_node->preg= 0;
						current_node->gest_at_doi= 0;
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
}

// Abortion: It returns an abortion indicater according to the gestation period (= "days_to_calv") of infected cattle. We assumed that there was no difference of abortion rate between TI and PI.
int _Abortion(int inf_stat, int days_to_calv) {
	int _Random_binomial(int n, double p);

	if (inf_stat == 1 || inf_stat == 4) {
		// Early pregnancy (abortion rate= 0.8)
		if (days_to_calv >= 238) {
			if (_Random_binomial(1, 0.8) == 1) {return 1;
			} else {return 0;}
		// Mid pregnancy (abortion rate= 0.25)
		} else if (days_to_calv >= 129 && days_to_calv < 238) {
			if (_Random_binomial(1, 0.25) == 1) {return 1;
			} else {return 0;}
		// Late pregnancy: No abortion
		} else {return 0;}
	} else {return 0;}
}

// Disease spread: It updates BVD infection status of each cattle.
void _Update_inf_stat(double *beta, int n_farms, int **group_size, animal_node *animal_list[]) {
	int _Random_binomial(int n, double p);
	double _Random_normal(double mu, double sd);
	void _Add_group_size(int **group_size, int farm_id, int demo_group, int inf_stat);
	void _Draw_group_size(int **group_size, int farm_id, int demo_group, int inf_stat);
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

		n_pi= group_size[i][4] + group_size[i][9];
		n_ti= group_size[i][1] + group_size[i][6];
			
		prev_pi= (double)n_pi / (double)n_herd;
		prev_ti= (double)n_ti / (double)n_herd;
			
		previous_node= animal_list[i];
		current_node= animal_list[i];
	
		// Probability of infection per demographic group
		p_inf= 1 - exp(-1.0 * ((*beta)*prev_pi + beta_ti*prev_ti));
		
		// Update transmission / recovery
		while (current_node != NULL) {
			// Recovered cattle
			if (current_node->demo_group == 99999 || current_node->demo_group == -99999 || current_node->inf_stat == 2) {
				previous_node= current_node;
				current_node= current_node->next_node;
			// Waning-off of Maternal Ab
			} else if (current_node->inf_stat == 3) {
				if (current_node->days_to_sr == 0) {
					_Draw_group_size(group_size, i, current_node->demo_group, current_node->inf_stat);
					current_node->inf_stat= 0;
					_Add_group_size(group_size, i, current_node->demo_group, current_node->inf_stat);
					previous_node= current_node;
					current_node= current_node->next_node;
				} else {
					current_node->days_to_sr -=1;
					previous_node= current_node;
					current_node= current_node->next_node;
				}
			// Transmission
			} else if (current_node->inf_stat == 0) {
				if (_Random_binomial(1, p_inf) == 1) {
					_Draw_group_size(group_size, i, current_node->demo_group, current_node->inf_stat);
					current_node->inf_stat= 1;
					current_node->days_to_sr= (int)_Random_normal(14.0, 2.0);
					if (current_node->preg == 1) {current_node->gest_at_doi= current_node->days_to_calv;}
					_Add_group_size(group_size, i, current_node->demo_group, current_node->inf_stat);
					previous_node= current_node;
					current_node= current_node->next_node;
				} else {
					previous_node= current_node;
					current_node= current_node->next_node;
				}
			// Recovery of TIs
			} else if (current_node->inf_stat == 1) {
				// Abortion of TI
				if (current_node->preg == 1 && _Abortion(current_node->inf_stat, current_node->days_to_calv) == 1) {
					current_node->preg= 0;
					current_node->days_to_calv= 0;
					current_node->gest_at_doi= 0;
				}
				// Recovery from TI (95% in 3 weeks): "gest_at_doi" SHOULDN'T BE 0 BY THIS STAGE
				if (current_node->days_to_sr == 0) {
					_Draw_group_size(group_size, i, current_node->demo_group, current_node->inf_stat);
					current_node->inf_stat= 2;
					_Add_group_size(group_size, i, current_node->demo_group, current_node->inf_stat);
					previous_node= current_node;
					current_node= current_node->next_node;
				// Not yet recovered
				} else {
					current_node->days_to_sr -=1;
					previous_node= current_node;
					current_node= current_node->next_node;
				}
			// Mortality of PIs
			} else {
				// PI survived
				if (_Random_binomial(1, 1.0 - lambda_pi) == 1) {
					// Abortion of PI
					if (current_node->preg == 1 && _Abortion(current_node->inf_stat, current_node->days_to_calv) == 1) {
						current_node->preg= 0;
						current_node->days_to_calv= 0;
						current_node->gest_at_doi= 0;
					}
					previous_node= current_node;
					current_node= current_node->next_node;
				// Death of PI
				} else {
					previous_node->next_node= current_node->next_node;
					_Draw_group_size(group_size, i, current_node->demo_group, current_node->inf_stat);
					free(current_node);
					current_node= previous_node->next_node;
				}
			}
		}
	}
}

// Sample the heifers for the test and run ELISA
void _Sample_serostatus(animal_node *animal_list[], int *pnt_disease[], int *n_sample) {
	animal_node *previous_node;
	animal_node *current_node;

	previous_node= animal_list[0];
	current_node= animal_list[0];

	int i= 0;
	while (current_node != NULL) {
		if (current_node->demo_group == 0) {
			if (i < (*n_sample)) {
				pnt_disease[i]= &current_node->inf_stat;
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
