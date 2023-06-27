// Prisoners Dilemma game on a small-world graph constructed from a square lattice 
// Some players are blocked to give their strategy (other players cannot adopt their strategy)

// standard include
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include <windows.h>
using namespace std;

// define priority classes
#define NORMAL_PRIORITY_CLASS       0x00000020
#define IDLE_PRIORITY_CLASS         0x00000040
#define HIGH_PRIORITY_CLASS         0x00000080
#define REALTIME_PRIORITY_CLASS     0x00000100

// define parameters
#define L           100     /* lattice size                   */
#define SIZE        (L*L)    /* number of sites                */
#define MC_STEPS    30000   /* run-time in MCS     */
//#define r               /* temptation to defect */
#define K           0.1      /* temperature */
#define Q           0.0      /* Q portion of links are rewired */
#define NAMEOUT     "K4b075r5Q2"
#define RANDOMIZE   3145215
#define delta       0.3


int C,D,loner,Bots;
double r;
int Time=10;
double rho; //rho is the fraction of bots


typedef int       tomb1[SIZE];
typedef long int  tomb3[SIZE][4];
typedef double    tomb4[SIZE];


tomb1 player_s;           /* matrix ,containing players strategies */
tomb3 player_n;           /* matrix, containing players neighbours */
tomb1 player_type;           /* matrix ,containing players type */



void prodgraph(void);      /* creates host graph                    */
void initial(void);        /* initial state                         */
void game(void);
void tongji(void);

ofstream outfile1;
ofstream outfile2;


//以下是随机数产生模块，不用管它,直接用就行，用randf()可以直接产生0-1满足均匀分布的随机数，randi(x),产生0---x-1的随机整数
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
static int mti=NN+1; /* mti==NN+1 means mt[NN] is not initialized */
void sgenrand(unsigned long seed)
{int i;
 for (i=0;i<NN;i++) {mt[i] = seed & 0xffff0000; seed = 69069 * seed + 1;
                     mt[i] |= (seed & 0xffff0000) >> 16; seed = 69069 * seed + 1;
  }
  mti = NN;
}
void lsgenrand(unsigned long seed_array[])
{ int i; for (i=0;i<NN;i++) mt[i] = seed_array[i]; mti=NN; }
double genrand() 
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    if (mti >= NN) 
    {
        int kk;
        if (mti == NN+1) sgenrand(4357); 
        for (kk=0;kk<NN-MM;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MM] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<NN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MM-NN)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[NN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[NN-1] = mt[MM-1] ^ (y >> 1) ^ mag01[y & 0x1];
        mti = 0;
    }  
    y = mt[mti++]; y ^= TEMPERING_SHIFT_U(y); y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C; y ^= TEMPERING_SHIFT_L(y);
    return y;  
}

double randf(){ return ( (double)genrand() * 2.3283064370807974e-10 ); }
long randi(unsigned long LIM){ return((unsigned long)genrand() % LIM); }

/********************** END of RNG ************************************/

//
//void initial(void)
//{
//	 int i,j;
//     TC=0;
//	 FC=0;
//	 D=0;
//	 ED=0;
//    for (i=0; i<L; i++)
//	{ 
//      for(j=0;j<L;j++){
//      	
//		player_s[i*L+j]=(int)randi(4);
//		
////		cout<<player_s[i*L+j]<<endl;
////	    if(0<=i && i<25) player_s[i*L+j]=0; 
////		if(25<=i && i<50) player_s[i*L+j]=1;
////		if(50<=i && i<70) player_s[i*L+j]=2;
////        if(75<=i && i<100) player_s[i*L+j]=3;
////        		
//		if(player_s[i*L+j]==0)  TC++;
//		else if(player_s[i*L+j]==1)  FC++;
//		else if(player_s[i*L+j]==2)  D++;
//		else ED++;
//	  
//	  }	
//	}
//}

void initial(void)
{
	 int i,j;
     C=0;
	 D=0;
	 loner=0;
	 Bots=0;
    for (i=0; i<SIZE; i++)
	{ 
	   // 
		if(randf()<rho){
		   player_type[i]=0; //bots
		}
		else player_type[i]=10; //normal people

		if(player_type[i]==0){
			player_s[i]=2;//the action of bots player
		}
		else player_s[i]=randi(3); // normal people have three actions: 0:C 1:D, and 2:loner 
		
		
		if(player_s[i]==0 && player_type[i]==10)  C++;
		else if(player_s[i]==1 && player_type[i]==10)  D++;
		else if(player_s[i]==2 && player_type[i]==10) loner++;
		else Bots++;

	  
	  }	

}

	
	   
// creates first a square grid graph and then rewires Q links 
void prodgraph(void)             
{
	int nneighbor, iu, ju, neighbor1, neighbor2;
	long int rewire, first, player1,player2,player3,MCe;
	double x;
	int i,j;
   

	// set up an initial square lattice
	for(i=0; i<L; i++)                     
	{
		for(j=0; j<L; j++)
		{ 
			// the first player
			player1 = L * j + i;             

			// and its four nearest neighbors
			iu = i + 1;         
			ju = j;     
			if (iu==L) iu = 0;
			player2 = L * ju + iu;  
			player_n[player1][0] = player2;

			iu = i;             
			ju = j + 1; 
			if (ju==L) ju = 0;
			player2 = L * ju + iu;  
			player_n[player1][1] = player2;

			iu = i - 1;         
			ju = j;     
			if (i==0) iu = L - 1;
			player2 = L * ju + iu;  
			player_n[player1][2] = player2;

			iu = i;             
			ju = j - 1; 
			if (j==0) ju = L - 1;
			player2 = L * ju + iu;  
			player_n[player1][3] = player2;
		}
	}

	// the rewiring starts - Q portion of joints is chosen randomly
	x = (double) (1.0 - Q);
	x = - log ( x );
	x = x*(2*SIZE);
	MCe = (long int) x;
	

	first = randi(SIZE);
	nneighbor = randi(4);
	player1 = player_n[first][nneighbor];
	neighbor1 = (nneighbor+2) % 4;

	for(rewire=0; rewire<MCe-1; rewire++)
	{
		do
		player2 = randi(SIZE);
		while(player2==first || player2==player1);
		do
		{
		neighbor2 = randi(4);
		player3 = player_n[player2][neighbor2];
		}
		while(player3==first || player3==player1);
		player_n[player1][neighbor1] = player2;
		player_n[player2][neighbor2] = player1;
		player1 = player3;
		if (player_n[player3][0]==player2)   neighbor1=0;
		if (player_n[player3][1]==player2)   neighbor1=1;
		if (player_n[player3][2]==player2)   neighbor1=2;
		if (player_n[player3][3]==player2)   neighbor1=3;
	 }	

	player_n[player1][neighbor1] = first;
	player_n[first][nneighbor] = player1;
  
	//cout<<player1<<'\t'<<player_n[player1][neighbor1]<<endl;
	
}

void game(void)
{
    int i,j,k;
	int strat1,strat2;
	int player1,player2;
	double P1,P2,F1,F2;
    int player11;
    int strat11;
    int suiji;                 //随机数
	double p,dP;
	double qqqq,w;

//	ofstream outfile3("细分.txt",ios::out|ios::app);
//
//		 if(!outfile3)
//		{
//        cout<<"不能打开此文件!";
//        exit(1);
//		}         
          

	for (i=0; i<SIZE; i++)
		{
			player1 = (int) randi(SIZE);     // randomly select a player1
		  
			strat1 = player_s[player1];   // the action of player1
			
			P1=0;
			
			F1=0;

			// the payoff of player1
			if(strat1==0) //if player1 chooses cooperation
			{
			   for(j=0;j<4;j++)
			   {
			    player11=player_n[player1][j];
				strat11=player_s[player11];
                  
		
				 if(strat11==0)
				  P1 +=1;
			
				 if(strat11==1)
				  P1 +=(-r);
				 
				 if(strat11==2)
				 P1 += delta;
				 
			   }
			}
			else if(strat1==1)  // if player1 chooses Defection
            {
			   for(j=0;j<4;j++)
			   {
			    player11=player_n[player1][j];
				strat11=player_s[player11];

				 if(strat11==0)
				  P1 +=(1+r);
				 
				 if(strat11==1)
				  P1 +=0;
				 
				 if(strat11==2)
				{
				 	P1 +=delta;			 
				 }

			   }
			}

			else if(strat1==2) //if player chooses loner strategy
			{

			  for(j=0;j<4;j++)
			   {
			    player11=player_n[player1][j];
				strat11=player_s[player11];

				 if(strat11==0)
				  P1 +=delta;
				 
				 if(strat11==1)
				  P1 +=delta;
				 if(strat11==2)
				   {
                  P1 +=delta;			 
				   }
				 
			   }
			}
		
			F1=P1;


		
	
//*****************************************************随机选择邻居
			suiji=randi(4);       
			


			player2 = player_n[player1][suiji];  //a randomly selected one neighbor:player2 

			strat2 = player_s[player2]; // the strategy of player2

			P2=0;
			
			F2=0;

			if(strat2==0) // if player2 chooses cooperation
			{
			   for(j=0;j<4;j++)
			   {

			    player11=player_n[player2][j];
				strat11=player_s[player11];
                  
		
				 if(strat11==0)
				 {
				  P2 +=1;
				 }
				 if(strat11==1)
				 {
				  P2 +=(-r);
				 }
				 if(strat11==2)
				 {
					 P2 +=delta;
				 }
			   }
			
			}
			else if(strat2==1) //if player2 chooses defection
            {
				for(j=0;j<4;j++)
			   {
			    player11=player_n[player2][j];
				strat11=player_s[player11];

				 if(strat11==0)
				 {
				  P2 +=(1+r);
				 }
				 if(strat11==1)
				 {
				  P2 +=0;
				 }
				 if(strat11==2)
				 {
				  P2 +=delta;
				 }
			   }
			}
			else if(strat2==2) // if player2 chooses loner
			{
				for(j=0;j<4;j++)
			   {
			    player11=player_n[player2][j];
				strat11=player_s[player11];

				 if(strat11==0)
				 {
				  P2 += delta;
				 }
				 if(strat11==1)
				 {
				  P2 += delta;
				 }
				 if(strat11==2)
				 {
				  P2 += delta;
				 }
			   }
			}
		
			F2=P2;

//***************************************************邻居
			
//          if(randf()<pow(10,-5)){ //mutation
          if(randf()<0){
          	while(1){
          	  int choice=randi(3);
          	  if(player_s[player1]!=choice&& player_type[player1]==10) 
              player_s[player1]=choice;
              break;
            //  cout<< player_s[player1]<<"\t"<<choice<<endl;
			  }

		   }
		else{
			
		   dP=F1-F2;
           p=1/(1+exp(dP*K)); 
           qqqq=randf();
           if(qqqq<p&&strat1!=strat2 && player_type[player1]==10)
		   player_s[player1]=strat2;
		   
	      } 



		 
//   static long mm = 1;
//	           if(mm%900==0) {
//		           cooperator=defector=punishment=0;
//
//		           for(int ii=0;ii<SIZE;ii++)
//						{
//	                        if(player_s[ii]==0)   cooperator++;
//				            if(player_s[ii]==1) defector++;
//                            if(player_s[ii]==2) punishment++;
//						}
//
//			      double m,xxx,YYY,ZZZ;
//		          m=(double)mm/SIZE;
//			      xxx=(double)cooperator/SIZE;
//				  YYY=(double)defector/SIZE;
//				  ZZZ=(double)punishment/SIZE;
//
//	    //cout<<m<<'\t'<<xxx<<'\t'<<endl;
//		          outfile3<<m<<'\t'<<xxx<<'\t'<<YYY<<'\t'<<ZZZ<<endl;
//			   } 
//	           mm += 1;
//	//if(mm==(M*M+1)) mm=1;  

         
		} //end of the SIZE  

}

           
void tongji(void)
{
  int i;
  C=0;
  D=0;
  loner=0;
  Bots=0;

  for(i=0;i<SIZE;i++)
  {
   if(player_s[i]==0 &&player_type[i]==10) C++;
   else if(player_s[i]==1&&player_type[i]==10) D++;
   else if(player_s[i]==2&&player_type[i]==10) loner++;
   else Bots++;

  }

}

void Adj(int t)
{
	
	ofstream outfile;
	char fname[20];
	sprintf(fname,"%d",t);
	strcat(fname,"adj.txt");
	outfile.open(fname);
	if(outfile.is_open())
	{
		for (int i=0;i<SIZE;i++)
		{
			for (int j=0;j<4;j++)
			{
				
				outfile<<i<<"\t"<<player_n[i][j]<<endl; 
			}
		}
		outfile.close(); 
	}
	else
	{
		cout<<"cannot open file !"<<endl;
	}
 }
 
 void Snapshot(int t)
{
	
	ofstream outfile;
	char fname[20];
	sprintf(fname,"%d",t);
	strcat(fname,".txt");
	outfile.open(fname);
	if(outfile.is_open())
	{
		for (int i=0;i<L;i++)
		{
			for(int j=0;j<L; j++){
			outfile<<player_s[i*L+j]+player_type[i*L+j]<<"\t"; 
			}
			outfile<<endl;				
		}
		outfile.close(); 
	}
	else
	{
		cout<<"cannot open file !"<<endl;
	}
 }
 
// the main program
int main()
{
	int steps,try_time;
    double aa,x,XX,x1,x2,x3,XX1,XX2,XX3,bb,cc,dd,XX_t,XX1_t,XX2_t,XX3_t;
    

	outfile1.open("frequency.txt");
	outfile2.open("average.txt");

    if(!outfile1)
	{
	 cout<<"can not open";
	 abort();
	}
   
	if(!outfile2)
	{
	 cout<<"can not open";
	 abort();
	}

	// initialize the random number generation
	sgenrand(RANDOMIZE);

	prodgraph();

//	Adj(steps==0);



	// begins the mc steps
		for(r=0.7;r<1.001;r=r+0.01) //dilemma strength
		  {  
		  	for(rho=0.46;rho<0.501;rho=rho+0.01){ //fraction of bots player
		  	
		  		XX_t=0;
				XX1_t=0;
				XX2_t=0;
				XX3_t=0;
			
		  	for(try_time=0;try_time<Time;try_time++){
		  	aa=0;
			bb=0;
			cc=0;
			dd=0;
			

			initial();
			for (steps=0; steps<MC_STEPS; steps++)
			{
			 x=(double)C/(SIZE-Bots);
		     x1=(double)D/(SIZE-Bots);
		     x2=(double)loner/(SIZE-Bots);;
		     x3=(double)Bots/SIZE;;
		     
			 game();
			 tongji();
//			 if(steps%1==0){
//			 				 outfile1<<steps<<'\t'<<x<<'\t'<<x1<<'\t'<<x2<<'\t'<<x3<<endl;
//			 cout<<steps<<'\t'<<x<<'\t'<<x1<<'\t'<<x2<<'\t'<<x3<<endl;
//			 }

		
	
//         if (steps==0||steps==5||steps==20||steps==50)
//	     {
//	       Snapshot(steps);
//	     }

		     if(steps>24999)
			   {
			   aa+=x;
			   XX=(double)aa/5000;
			   }
		      if(steps>24999)
			   {	
	           bb+=x1;
			   XX1=(double)bb/5000;	
			   }
		     if(steps>24999)
			   {
	     	   cc+=x2;
			   XX2=(double)cc/5000;
			   }
//			 if(steps>24999)
//			   {
//	     	   dd+=x3;
//			   XX3=(double)dd/5000;
//			   }
            if(steps>=2000){
            	if(x==1){
			   XX=1;
			   XX1=0;
			   XX2=0;
			   XX3=x3;
			   break;
			   }
			   if(x1==0){
			   XX=x;
			   XX1=0;
			   XX2=x2;
			   XX3=x3;
			   break;
			   }
			   if(x2==1){
			   XX=0;
			   XX1=0;
			   XX2=1;
			   XX3=x3;
			   break;
			   }
////			   if(x3==1){
////			   XX=0;
////			   XX1=0;
////			   XX2=0;
////			   XX3=1;
////			   break;
////			   }
			}			   
			}//end of the MC steps
			outfile1<<try_time<<'\t'<<r<<"\t"<<rho<<"\t"<<x<<'\t'<<x1<<'\t'<<x2<<'\t'<<x3<<endl;
			 cout<<try_time<<'\t'<<x<<'\t'<<x1<<'\t'<<x2<<'\t'<<x3<<endl;
			 
			 XX_t+=XX;
			 XX1_t+=XX1;
			 XX2_t+=XX2;
			 XX3_t+=XX3;
			  }//end oftry
		  		   
			   
			    outfile2<<r<<'\t'<<rho<<'\t'<<XX_t/Time <<'\t'<<XX1_t/Time<<'\t'<<XX2_t/Time<<'\t'<<XX3_t/Time<<endl;
			    cout<<r<<'\t'<<rho<<'\t'<<XX_t/Time <<'\t'<<XX1_t/Time<<'\t'<<XX2_t/Time<<'\t'<<XX3_t/Time<<endl;
		  	
		}//end of rho 
	} // end of b
	
  
	outfile1.close();
	outfile2.close();
	
	return 0;
}	
