#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>

using namespace std;
#define L 100
#define N L*L
#define RANDOMIZE   3145215
#define str_num  2
#define neig_num  8

int neighbors[N][neig_num], strategy[N],punishment_strategy[N];
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


void find_neig(void)
{
	for(int i = 0 ; i<N ; i++)
	{
		neighbors[i][0]=i-L; //up 
		neighbors[i][1]=i+L; //down
		neighbors[i][2]=i-1; //left
		neighbors[i][3]=i+1; //right
		neighbors[i][4]=i-L-1; //upper left 
		neighbors[i][5]=i-L+1; //upper right
		neighbors[i][6]=i+L-1; //down left
		neighbors[i][7]=i+L+1; //down right
		if (i<L)
		{
			neighbors[i][0]= i + L * (L-1);
			neighbors[i][4]= i + L * (L-1)-1;
			neighbors[i][5]= i + L * (L-1)+1;
		}
		if (i > L * (L - 1) -1)
		{
			neighbors[i][1]= i - L * (L-1);
			neighbors[i][6]= i - L * (L-1)-1;
			neighbors[i][7]= i - L * (L-1)+1;
		}
		  
		if (i%L == 0)      
		{
			neighbors[i][2]= i + L - 1 ;
			neighbors[i][4]=i-1;
			neighbors[i][6]=i+2*L-1; 
		}     
		if (i%L == L - 1)     
		{
			neighbors[i][3]= i - L + 1 ; 
			neighbors[i][5]= i-2*L+1; 
			neighbors[i][7]= i + 1;
		}   
		if (i == 0)                  neighbors[i][4]= L*L-1; 
		else if (i == L-1)           neighbors[i][5]= L*(L-1);  
		else if (i == L*(L-1))       neighbors[i][6]= L-1;
		else if (i == L*L-1)         neighbors[i][7]= 0;  
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
	for (int i = 0; i < neig_num; i++) {
		neig = neighbors[x][i];
		pay += payoff_matrix[strategy[x]][strategy[neig]];
		if (punishment_strategy[x] == 0)//If the player chooses to punish
		{
			if (strategy[neig] == 1)
			{
				a=a+1;
			}
		}
		if (strategy[x] == 1)//If the player at the location chooses to betray
		{
			if (punishment_strategy[neig] == 0)
			{
				b=b+1;
			}
		}
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

#define loop 1
//double record_loop[loop][4];

//double mean_temp_cp,mean_temp_cn,mean_temp_dp,mean_temp_dn;
//void cal_dev(void){
	//mean_temp_cp=0;
	//mean_temp_cn=0;
	//mean_temp_dp=0;
//	mean_temp_dn=0;
	//for(int i=0; i<loop; i++) {
//		mean_temp_cp += record_loop[i][0];
//		mean_temp_cn += record_loop[i][1];
//		mean_temp_dp += record_loop[i][2];
//		mean_temp_dn += record_loop[i][3];
//	}
//	mean_temp_cp /= (double) loop;
//	mean_temp_cn /= (double) loop;
//	mean_temp_dp /= (double) loop;
//	mean_temp_dn /= (double) loop;
//}




int main(void) {
    int Round = 100000;
    //int mont_step = N ;
    sgenrand(time(0)); // Sets the current time as a random seed, affects randi() and randf().

    printf("*****start*****\n");

    FILE *Fc = fopen("cost and fine R=0.2 delta=0.03 beta=0.3 direct_payoff_time.csv", "w");  // 打开文件以保存合作率数据
    R = 0.2 ;
    delta=0.03;
    beta=0.3;
    init_game(R);
    for(int j = 0; j < Round; j++){
    	round_game();// play game
        cal_data();//caculate fractions of strategies
        double a = 0, b = 0;
        double c = 0, d = 0;
	double payoff_cost = 0;
	double payoff_fine = 0;
	double cost = 0;
	double fine = 0;
	for (int y = 0; y < N; y++) {
		
		for (int i = 0; i < neig_num; i++) {
			int neig;
			neig = neighbors[y][i];
			if (strategy[y]==0&&punishment_strategy[y] == 0&&strategy[neig]==1){
			    a=a+1;
			    for (int i = 0; i < neig_num; i++) {
			    	int neig2;
			    	neig2 = neighbors[y][i];
			    	if (strategy[neig2] == 1)b=b+1;//D
			    }  		
		    }
		   
		    else if (strategy[neig]==0&&punishment_strategy[neig]==0&&strategy[y]==1){
		    	c=c+1;
		    	for (int i = 0; i < neig_num; i++) {
				    int neig2;
				    neig2 = neighbors[y][i];
			  	    if (strategy[neig2]==0&&punishment_strategy[neig2] == 0)d=d+1;//CP
			    }
		    }
		}
	}
		

    payoff_cost = (b*delta)/a;
	payoff_fine = (d*beta)/c;
	cost = b/a;
	fine = d/c;

        printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",j,R,delta,beta,data_out[0],data_out[1],data_out[2],data_out[3],payoff_cost,payoff_fine,cost,fine);
        
         fprintf(Fc,"%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",j,R,delta,beta,data_out[0],data_out[1],data_out[2],data_out[3],payoff_cost,payoff_fine,cost,fine);

	}
    
    fclose(Fc);
    printf("*****done*****");
    return 0;
}
