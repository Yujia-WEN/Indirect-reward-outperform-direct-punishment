#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstring>

using namespace std;
#define L 100
#define N L*L
#define RANDOMIZE   3145215
#define str_num  2
#define neig_num  8

int neighbors[N][neig_num], strategy[N],punishment_strategy[N];
int second_neig[N][16] = {0};

double payoff_matrix[str_num][str_num];

double K = 0.1;
double R = 0.2;
double delta;
double beta;

//The following is the random number generation module, use randf() to directly generate 0-1 random numbers that satisfy the uniform distribution, randi(x), generate 0---x-1 random integers
/*************************** RNG procedures ****************************************/
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[NN]; /* the array for the state vector  */
static int mti = NN + 1; /* mti==NN+1 means mt[NN] is not initialized */
void sgenrand(unsigned long seed) {
	int i;
	for (i = 0; i < NN; i++) {
		mt[i] = seed & 0xffff0000;
		seed = 69069 * seed + 1;
		mt[i] |= (seed & 0xffff0000) >> 16;
		seed = 69069 * seed + 1;
	}
	mti = NN;
}
void lsgenrand(unsigned long seed_array[]) {
	int i;
	for (i = 0; i < NN; i++) mt[i] = seed_array[i];
	mti = NN;
}
double genrand() {
	unsigned long y;
	static unsigned long mag01[2] = { 0x0, MATRIX_A };
	if (mti >= NN) {
		int kk;
		if (mti == NN + 1) sgenrand(4357);
		for (kk = 0; kk < NN - MM; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + MM] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		for (; kk < NN - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (MM - NN)] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		y = (mt[NN - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[NN - 1] = mt[MM - 1] ^ (y >> 1) ^ mag01[y & 0x1];
		mti = 0;
	}
	y = mt[mti++];
	y ^= TEMPERING_SHIFT_U(y);
	y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
	y ^= TEMPERING_SHIFT_L(y);
	return y;
}

double randf() {
	return ((double)genrand() * 2.3283064370807974e-10);
}
long randi(unsigned long LIM) {
	return((unsigned long)genrand() % LIM);
}

/********************** END of RNG ************************************/

int temp_net[L+4][L+4];
void find_neig(void)
{
	memset(temp_net, 0, sizeof(temp_net));
	memset(neighbors, 0, sizeof(neighbors));
	memset(second_neig, 0, sizeof(second_neig));
	
	int intera = 2;
	int count=0;
	for(int i=intera; i<L+intera; i++){
		for(int j=intera; j<L+intera; j++){
			temp_net[i][j] = count;
			count ++;
		}
	}
	
	int temp_x;
	for(int i=0; i<intera; i++){
		for(int j=intera; j<L+intera; j++){
			temp_net[i][j] = temp_net[i+L][j];
			temp_net[i+L+intera][j] = temp_net[intera+i][j];
			temp_net[j][i] = temp_net[j][i+L];
			temp_net[j][i+L+intera] = temp_net[j][intera+i];
		}
	}
	for(int i=0; i<intera; i++){
		for(int j=0; j<intera; j++){
			temp_net[i][j]  = temp_net[i][j+L];
			temp_net[i][j+L+intera]  = temp_net[i][j+intera];
			temp_net[i+L+intera][j]  = temp_net[i+L+intera][j+L];
			temp_net[i+L+intera][j+L+intera]  = temp_net[i+L+intera][j+intera];
		}
	}
	
	
	intera = 1; 
	int center;
	int temp_node;
	int neig_count;
	for(int i=intera; i<L+intera; i++){
		for(int j=intera; j<L+intera; j++){
			center = temp_net[i][j];
			neig_count = 0;
			for(int new_i=i-intera; new_i<i+intera+1; new_i++){
				for(int new_j=j-intera; new_j<j+intera+1; new_j++){
					temp_node = temp_net[new_i][new_j];
					if(temp_node != center){
						neighbors[center][neig_count] = temp_node;
						neig_count ++;
					}
				}
			}
		}
	}
	
	
	intera = 2;
	for(int i=intera; i<L+intera; i++){
		for(int j=intera; j<L+intera; j++){
			center = temp_net[i][j];
			neig_count = 0;
			
			for(int new_i=i-intera; new_i < i+intera+1; new_i++){
				if( new_i == i-intera or new_i == i+intera ){
					for(int new_j=j-intera; new_j<j+intera+1; new_j++){
						
						temp_node = temp_net[new_i][new_j];
						second_neig[center][neig_count] = temp_node;
						neig_count ++;
					}
					
				}
				
				else{
					temp_node = temp_net[new_i][j-intera];
					second_neig[center][neig_count] = temp_node;
					neig_count ++;
					
					temp_node = temp_net[new_i][j+intera];
					second_neig[center][neig_count] = temp_node;
					neig_count ++;
				} 
			}
		}
	}
}

//Initialize the payoff matrix
void init_game(double R)
{
	find_neig();
	payoff_matrix[0][0] = 1;
	payoff_matrix[0][1] = (double)(-R);
	payoff_matrix[1][0] = (double) (1 + R);
	payoff_matrix[1][1] = 0;


//Initialize the strategy and punishment strategy for the player
	for (int i = 0; i < N; i++) {
		strategy[i] = randi(str_num);
		punishment_strategy[i] = randi(str_num);
	}
}

// Calculate payoff
double cal_payoff(int x) {
	int a = 0;
	int b = 0;
	double pay = 0;
	int neig;
	int neigneig;
	for (int i = 0; i < neig_num; i++) {
		neig = neighbors[x][i];
		pay += payoff_matrix[strategy[x]][strategy[neig]];

	}
	for(int j = 0; j < 16; j++) {
		neigneig=second_neig[x][j];
		if(strategy[neigneig]==1&&punishment_strategy[x]==0)a++;
		if(strategy[x]==1&&punishment_strategy[neigneig]==0)b++;
	}
		
		
    pay = pay-a*delta-b*beta;
    return pay;	
}

void learn_strategy(int x) {
	double epsil = 1e-5;
	if (randf() < epsil)
	{
		strategy[x] = randi(str_num);
		punishment_strategy[x] = randi(str_num);
	}
	else{
		int neig = neighbors[x][randi(neig_num)];
		double x_r = cal_payoff(x);
	    double n_r = cal_payoff(neig);

	    double prob = (double)1 / (1 + exp((x_r - n_r) / K));
	/*Femi Function  K=0.1 means the selection intensity, the smaller the selection intensity, the greater the selection intensity.
	 K=0.1 means that the learning probability is between 70% and 80%. The reason why there is no 100% is because human society is not completely rational */
	//prob imitates the neighbor's strategy, 1-prob maintains the original strategy
	    double p = randf();
	    if (p < prob) 
	    {
	    	strategy[x] = strategy[neig];
		    punishment_strategy[x] = punishment_strategy[neig];
		}
	}

	
}

void round_game(void) {
	int center;
	//Asynchronous Simulation: Synchronous Simulation
	for (int i = 0; i < N; i++) {

		center = randi(N);
		learn_strategy(center);
	}
}

double data_out[4];// Create an empty array of length 4 to store the data
void cal_data() {
	int x = 0, y = 0, z = 0, w = 0;
	int nn = N;
	for (int i = 0; i < N; i++) {
		if (strategy[i] == 0 && punishment_strategy[i] == 0) x++;
	    else if (strategy[i] == 0 && punishment_strategy[i] == 1) y++;
	    else if (strategy[i] == 1 && punishment_strategy[i] == 0) z++;
	    else if (strategy[i] == 1 && punishment_strategy[i] == 1) w++;
	}
	data_out[0] = (double)x / nn;
	data_out[1] = (double)y / nn;
	data_out[2] = (double)z / nn;
	data_out[3] = (double)w / nn;// output cooperation rate
}

#define loop 10
double record_loop[loop][4];

double mean_temp_cp,mean_temp_cn,mean_temp_dp,mean_temp_dn;
void cal_dev(void){
	mean_temp_cp=0;
	mean_temp_cn=0;
	mean_temp_dp=0;
	mean_temp_dn=0;
	for(int i=0; i<loop; i++) {
		mean_temp_cp += record_loop[i][0];
		mean_temp_cn += record_loop[i][1];
		mean_temp_dp += record_loop[i][2];
		mean_temp_dn += record_loop[i][3];
	}
	mean_temp_cp /= (double) loop;
	mean_temp_cn /= (double) loop;
	mean_temp_dp /= (double) loop;
	mean_temp_dn /= (double) loop;
}





int main(void) {
	int Round = 30000;
	//int mont_step = N ;
	sgenrand(time(0));// Sets the current time as a random seed, affects randi() and randf().



	printf("*****start*****\n");


	FILE *Fc = fopen("r=0.2 deta=0.20 indirect.csv", "w"); //"w": mode, open for writing
    
	 for  (delta = 0.00; delta < 1.01; delta += 0.01){
    	   beta = 0.03;
    		                          
    		for(int lo =0; lo<loop; lo++){
    		    	                     
						init_game(R);
    		    double cp_temp = 0, dp_temp = 0;
    		    double cn_temp = 0, dn_temp = 0;
    		    
    			for (int i = 0; i < Round; i++) { // round

			    round_game();// play game
			    
			    cal_data();//caculate fractions of strategies
			    if (i>= Round - 2000){
			    	cp_temp = cp_temp + data_out[0];
			    	cn_temp = cn_temp + data_out[1];
			    	dp_temp = dp_temp + data_out[2];
			    	dn_temp = dn_temp + data_out[3];
				}  
				
				
				}
				record_loop[lo][0] = (double)cp_temp/2000;
				record_loop[lo][1] = (double)cn_temp/2000;
				record_loop[lo][2] = (double)dp_temp/2000;
				record_loop[lo][3] = (double)dn_temp/2000;
			
			
			printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n",R,delta,beta,(double)cp_temp/2000,(double)cn_temp/2000,(double)dp_temp/2000,(double)dn_temp/2000);
		}
		
		cal_dev();
		fprintf(Fc,"%f,%f,%f,%f,%f,%f,%f\n",R,delta,beta,mean_temp_cp,mean_temp_cn,mean_temp_dp,mean_temp_dn);
		
		
	}


    fclose(Fc);
	printf("*****done*****");
	return 0;
}
