#include<algorithm>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cstring>
#include<vector>
#include<ctime>
#include<cmath>
using namespace std;


#define Max_user    500000

#define Max_channel    10
#define Max_level      10

// arguments

const int phy_r[4] = {27000, 18000, 9000, 3000}; // physical rates(kbps)
const int s_type[4] = {4500, 2000, 1000, 300}; // service types 0:video, 1:FTP, 2:web, 3:VoIP


// utility arguments

double a1 = 50, b1 = 8;
double a2 = 5;
double a3 = 4.5;
double a4 = 140, b4=12;


int u_type[4]; // users' service type
double ml[4];
double sl[4];

int n_phy_g[4];

const double FT = 2; // frame time 0.4s
const double ST = 0.01; // slot time 0.01s
const int z = 75; // 15 slots each frame
int M;            // number of users
int N;            // number of request channels
int N_;           // number of got channels

double gama = 0.2;
double beta = 0.2;
double alpha;     

int s_t=0;//scheduling times


// math functions
int ceiling (double x) {
    return x - (int)x > 0 ? (int)x + 1 : (int)x;
}

int i_floor (double x) {
    return (int)x;
}


// get utility functions

double u1(double ur){//video
	
	return 1.0 / (1.0 + a1* exp(-1 * b1 * (ur / s_type[0] )));
}
double u2(double ur){//FTP
	return 1.0 - exp(-1 * a2 * (ur / s_type[1] ));
}

double u3(double ur){//web
	return 1.0 - exp(-1 * a3 * (ur / s_type[2] ));
}
double u4(double ur){//VoIP
	return 1.0 / (1.0 + a4* exp(-1 * b4 * (ur / s_type[3] )));
}

double (*uti[4])(double ur);


class flow_current {
public:
    
    int type;
    double got_r;//現在拿到多少rate 
    double max_r;
    double min_r; 
    int AP;
    int sector;
    long ID;
    int allocated_slot[Max_channel][z];
    int allo_time[z];//某個時間點的slot 
    double priority;
    double queue_data;
    double utility_c;//current utility 
    

};


class flow_all {
public:
    int ser_type; // service types 0:video, 1:FTP, 2:web, 3:VoIP
    long MU_ID;
    double rate; 
    double utility; 
    int sch_times;

};

int main () 
{
	uti[0] = u1;
	uti[1] = u2;
	uti[2] = u3;
	uti[3] = u4;

	int current_AP=0;
	int current_sector=0;
	vector<flow_all> f_all (Max_user);
	int num_MU=0;
	int type_sum=0;
	int SU_no_ch=0;
	double ave_u=0.0;
	int no_gch_f[4]={0,0,0,0};//number of blocking user
	for(int i=0 ; i<Max_user ; i++){ 
           f_all[i].MU_ID = -1;//if no this flow => -1
           f_all[i].rate = 0;
           f_all[i].utility = 0;
           f_all[i].sch_times = 0;
           f_all[i].ser_type = -1;
           
           
    } 
	

	// input
	printf("Please enter input file:\n");
	string filename;
    cin >> filename;
    ifstream fin(filename.c_str());
    ofstream fout("yu_algo_result_offpeak.txt");
    ofstream fout2("yu_algo_result_all_flow_offpeak.txt");
    
    while (fin >> current_AP) {
        fin >> current_sector;
        
        fin >> num_MU;
        
        cout << "AP " << current_AP << " sector " << current_sector << endl;
        fout << "AP " << current_AP << " sector " << current_sector << endl;
        vector<flow_current> f_cur (num_MU);  
        
        
        
        
        double t_r = 0;
        for(int i=0 ,j=0; i<4 ;i++){
                fin >> u_type[i];
                t_r += u_type[i]*s_type[i]*FT;
                
                for(int k=0 ; k< u_type[i] ; j++,k++){
                        fin >> f_cur[j].ID;
                        f_cur[j].type = i;//0~3
                        f_cur[j].AP = current_AP;
                        f_cur[j].sector = current_sector;
                        f_cur[j].queue_data = s_type[i]*FT;
                        
                }        
        }
        
        double t_g = 0;
            
        fin >> N_;
        
        int slot_result[N_][z];
        int slot_type[N_][z];
        
        
       
        for (int i = 0,k=0 ; i < 4 ; i++) {             
        	    fin >> n_phy_g[i];    	        
        	    t_g += n_phy_g[i] * phy_r[i] * z * ST;
        	    for(int j=0;j<n_phy_g[i];j++){
                        
                      for(int a=0;a<z;a++){                
                           slot_type[k][a]=i;
                           slot_result[k][a]=-1;//not allocated to any user             
                      }
                      k++;                               
                }
       	}
       	
      
        if(num_MU==0){
              //cout << "No any flow\n" << endl;
              //fout << "No any flow\n" << endl;
              continue;
        }
        s_t++;
        
        if(N_==0 && num_MU!=0 ){//no gch flow not insert to flow structure
              for(int i=0;i<4;i++)
                   no_gch_f[i]+=u_type[i];
              SU_no_ch++;
              cout << "No granted channel\n" << endl;
              fout << "No granted channel\n" <<  endl;
              continue;
        }
        
        if(t_r<t_g)
                t_r=t_g;
                	
		alpha = (t_g / t_r) * gama;
		
		for(int i=0 ;i<4; i++){
                if(beta==0)
                    ml[i]=50000;
                else
                    ml[i] = s_type[i] / beta;
                sl[i] = s_type[i] * alpha; 
        }
        for(int i=0 ; i<num_MU ; i++){ 
                f_cur[i].max_r = ml[f_cur[i].type];
                f_cur[i].min_r = sl[f_cur[i].type];
                f_cur[i].got_r = 0;
               
                f_cur[i].priority = 0;
                f_cur[i].utility_c = 0;
                f_all[f_cur[i].ID].MU_ID = f_cur[i].ID;
                f_all[f_cur[i].ID].ser_type = f_cur[i].type;
                
                f_all[f_cur[i].ID].sch_times ++; 
                
                for(int j=0;j<z;j++)
                        f_cur[i].allo_time[j]=0;
                for(int m=0;m<Max_channel;m++){
                    for(int n=0;n<z;n++){
                        f_cur[i].allocated_slot[m][n] =0; 
                    }
                }      
        } 
        
        
        
		//Phase1
		
		bool ini_ok = true;
        //for (int i = 0 ; i < num_MU ; ++i) { 
            
        for (int i = num_MU-1 ; i >=0 ; i--) {   //low rate flow first   
                bool found = false;
               
                for(int m=0;m<N_;m++){// high rate slot first
                //for(int m=N_-1;m>=0;m--){        
                    for(int n=0;n<z;n++){// if allocated   slot_type[m][n]=-1
                        if(slot_type[m][n]==-1 || f_cur[i].allo_time[n]==1)
                              continue;  
                              
                        //got min guarantee
                        if(f_cur[i].got_r >= f_cur[i].min_r && f_cur[i].got_r <= f_cur[i].max_r){
                              found=true;
                              break;           
                        }
                                   
                        if(f_cur[i].got_r +(( phy_r[slot_type[m][n]] *ST)/FT)<= f_cur[i].max_r){  	
                              f_cur[i].allocated_slot[m][n]=1;
                              f_cur[i].got_r+=( phy_r[slot_type[m][n]] *ST)/FT;
                              f_cur[i].queue_data -= phy_r[slot_type[m][n]] *ST;
                              slot_result[m][n]=f_cur[i].ID;
                              slot_type[m][n]=-1;
                              f_cur[i].allo_time[n]=1;
                        }
                        
                       
                        
                            
                    }   
                    if(found)
                        break;                        
               	}
                        
                if(!found){//someone constraint failed
                    ini_ok= false;
                    break;
                }			
        }
        
        	
        if (!ini_ok) {
				cout << "Phase1 failed" << endl;
				fout << "Phase1 failed" << endl;
				continue;
		}else{
                for (int i = 0 ; i < num_MU ; ++i){ 
                    f_all[f_cur[i].ID].rate += f_cur[i].got_r;
                    
                }
                    
        }
		
        cout << "Phase1 end" << endl;
		//Phase1 end
		
		//Phase 2
		//for(int m=N_-1;m>=0;m--){
		for(int m=0;m<N_;m++){// high rate slot first
               for(int n=0;n<z;n++){
                    if(slot_type[m][n]==-1 ){
                         continue;
                    }
                    double max_gradient = -1; 
                    double gradient=0;
                    int choose_flow =-1;
                       
                    for (int i = 0 ; i < num_MU ; ++i) { 
                         double r_promoted = f_cur[i].got_r + (( phy_r[slot_type[m][n]] *ST)/FT);
                         //if no queue data or max guarantee failed or time already allocated
                         if(f_cur[i].queue_data<= 0 || r_promoted > f_cur[i].max_r || f_cur[i].allo_time[n]==1){
                             continue;
                         }
                         gradient= (uti[f_cur[i].type](r_promoted) - uti[f_cur[i].type](f_cur[i].got_r ))/ (r_promoted- f_cur[i].got_r);
                         if(f_all[f_cur[i].ID].sch_times!=0 && f_all[f_cur[i].ID].utility!=0 ) 
                              f_cur[i].priority = 1/(f_all[f_cur[i].ID].utility/f_all[f_cur[i].ID].sch_times);
                         else
                              f_cur[i].priority = 100000;//max priority 
                         
                         
                         if(gradient >= max_gradient){
                              if(choose_flow==-1){//if haven't chosen any flow
                                   max_gradient = gradient;
                                   choose_flow = i;
                              }else{
                              
                                   if(gradient == max_gradient ){
                                        if(f_cur[i].priority > f_cur[choose_flow].priority){
                                             max_gradient = gradient;
                                             choose_flow = i;
                                             
                                        }
                                   }else{
                                         max_gradient = gradient;
                                         choose_flow = i;
                                   }         
                                            
                              }
                              
                         }
                         
                    }//for num_MU
                    if(choose_flow!=-1){
                         f_cur[choose_flow].allocated_slot[m][n]=1;
                         f_cur[choose_flow].got_r+=( phy_r[slot_type[m][n]] *ST)/FT;
                         f_all[f_cur[choose_flow].ID].rate +=( phy_r[slot_type[m][n]] *ST)/FT;
                         f_cur[choose_flow].queue_data -= phy_r[slot_type[m][n]] *ST;
                         slot_result[m][n]=f_cur[choose_flow].ID;
                         slot_type[m][n]=-1;
                         f_cur[choose_flow].allo_time[n]=1;
                    }else{
                         //cout<<"No one can use this slot ["<< m<<"]["<<n<<"]"<<endl;
                    }            
               }
        }
        double total_u=0;
        for (int i = 0 ; i < num_MU ; ++i) { 
               f_cur[i].utility_c = uti[f_cur[i].type](f_cur[i].got_r);
               f_all[f_cur[i].ID].utility += f_cur[i].utility_c;
              
               total_u+=f_cur[i].utility_c;
               cout<<"flow :" << f_cur[i].ID <<" type :"<<f_cur[i].type <<" uti :"<<f_cur[i].utility_c<< " rate :"<< f_cur[i].got_r <<endl;
               fout<<"flow :" << f_cur[i].ID <<" type :"<<f_cur[i].type <<" uti :"<<f_cur[i].utility_c<< " rate :"<< f_cur[i].got_r <<endl;
        }
        
        for(int m=0;m<N_;m++){
               for(int n=0;n<z;n++){
                    
                    fout<<"["<<slot_result[m][n]<<"] ";
               }
              
               fout<<endl;
        }
        cout<<"ave utility :" << total_u/num_MU <<endl;
        fout<<"ave utility :" << total_u/num_MU <<"\n"<<endl;
		cout << "Phase2 end" <<"\n"<< endl;
        
	}//while(AP)
	
	//output scheduling result
	double ave_rate[4]={0,0,0,0},ave_utility[4]={0,0,0,0};
	int num_count[4]={0,0,0,0};
	double rate_level[4][Max_level+1];
	double utility_level[4][Max_level+1];
	int ser_level_count[4][Max_level+1];
	int uti_level_count[4][Max_level+1];
	
	
	for(int i=0 ; i<4 ; i++){
            for(int j=0 ; j<=Max_level ;j++){
                    rate_level[i][j]=j*(s_type[i]/Max_level);
                    ser_level_count[i][j]=0;
                    utility_level[i][j]=(double)(1.0/Max_level)*j;
                    uti_level_count[i][j]=0;
                    
            }
    }
	
	for(int i=0 ; i<Max_user  ; i++){         
           if(f_all[i].ser_type!=-1){   
               cout<<"ID: " <<f_all[i].MU_ID<<" type: "<<f_all[i].ser_type<<" rate: " << f_all[i].rate<<" uti: "<< f_all[i].utility <<" sch times: "<<f_all[i].sch_times <<endl;
               fout2<<"ID: " <<f_all[i].MU_ID<<" type: "<<f_all[i].ser_type<<" rate: " << f_all[i].rate<<" uti: "<< f_all[i].utility <<" sch times: "<<f_all[i].sch_times <<endl;
               if(f_all[i].ser_type==0){//Video
                   ave_rate[0]+= (f_all[i].rate/f_all[i].sch_times);
                   ave_utility[0]+= (f_all[i].utility/f_all[i].sch_times);
                   num_count[0]++;
                   
                   for (int j=0;j<=Max_level ;j++){
                       
                       
                       if((f_all[i].rate/f_all[i].sch_times) > rate_level[0][Max_level-1]+((rate_level[0][Max_level]-rate_level[0][Max_level-1])/2)){
                            ser_level_count[0][Max_level]++;
                            
                            break;
                       }
                       
                       if((f_all[i].rate/f_all[i].sch_times) <=rate_level[0][j]+((rate_level[0][j+1]-rate_level[0][j])/2)){
                                                            
                            ser_level_count[0][j]++;
                            
                            break;
                       }
                      
                       
                   }
                   for (int j=0;j<=Max_level ;j++){
                       
                       if((f_all[i].utility/f_all[i].sch_times) > utility_level[0][Max_level-1]+((utility_level[0][Max_level]-utility_level[0][Max_level-1])/2)){
                            uti_level_count[0][Max_level]++;
                            
                            break;
                       }
                       
                       if((f_all[i].utility/f_all[i].sch_times) <=utility_level[0][j]+((utility_level[0][j+1]-utility_level[0][j])/2)){
                            uti_level_count[0][j]++;
                            
                            break;
                       }
                       
                   }
               }else if(f_all[i].ser_type==1){//FTP
                   ave_rate[1]+= (f_all[i].rate/f_all[i].sch_times);
                   ave_utility[1]+= (f_all[i].utility/f_all[i].sch_times);
                   num_count[1]++;
                   
                   for (int j=0;j<=Max_level ;j++){
                       if((f_all[i].rate/f_all[i].sch_times) > rate_level[1][Max_level-1]+((rate_level[1][Max_level]-rate_level[1][Max_level-1])/2)){
                            ser_level_count[1][Max_level]++;
                            
                            break;
                       }
                       
                       if((f_all[i].rate/f_all[i].sch_times) <=rate_level[1][j]+((rate_level[1][j+1]-rate_level[1][j])/2)){
                            ser_level_count[1][j]++;
                            
                            break;
                       }
                       
                   }
                   for (int j=0;j<=Max_level ;j++){
                       
                       if((f_all[i].utility/f_all[i].sch_times) > utility_level[1][Max_level-1]+((utility_level[1][Max_level]-utility_level[1][Max_level-1])/2)){
                            uti_level_count[1][Max_level]++;
                            
                            break;
                       }
                       
                       if((f_all[i].utility/f_all[i].sch_times) <=utility_level[1][j]+((utility_level[1][j+1]-utility_level[1][j])/2)){
                            uti_level_count[1][j]++;
                            
                            break;
                       }
                       
                   }
               }else if(f_all[i].ser_type==2){//Web
                   ave_rate[2]+= (f_all[i].rate/f_all[i].sch_times);
                   ave_utility[2]+= (f_all[i].utility/f_all[i].sch_times);
                   num_count[2]++;
                   
                   for (int j=0;j<=Max_level ;j++){
                       if((f_all[i].rate/f_all[i].sch_times) > rate_level[2][Max_level-1]+((rate_level[2][Max_level]-rate_level[2][Max_level-1])/2)){
                            ser_level_count[2][Max_level]++;
                            
                            break;
                       }
                       
                       if((f_all[i].rate/f_all[i].sch_times) <=rate_level[2][j]+((rate_level[2][j+1]-rate_level[2][j])/2)){
                            ser_level_count[2][j]++;
                           
                            break;
                       }
                       
                   }
                   for (int j=0;j<=Max_level ;j++){
                       
                       if((f_all[i].utility/f_all[i].sch_times) > utility_level[2][Max_level-1]+((utility_level[2][Max_level]-utility_level[2][Max_level-1])/2)){
                            uti_level_count[2][Max_level]++;
                            
                            break;
                       }
                       
                       if((f_all[i].utility/f_all[i].sch_times) <=utility_level[2][j]+((utility_level[2][j+1]-utility_level[2][j])/2)){
                            uti_level_count[2][j]++;
                            
                            break;
                       }
                       
                   }
               }else if(f_all[i].ser_type==3){//VoIP
                   ave_rate[3]+= (f_all[i].rate/f_all[i].sch_times);
                   ave_utility[3]+= (f_all[i].utility/f_all[i].sch_times);
                   num_count[3]++;
                   
                   for (int j=0;j<=Max_level ;j++){
                       if((f_all[i].rate/f_all[i].sch_times) > rate_level[3][Max_level-1]+((rate_level[3][Max_level]-rate_level[3][Max_level-1])/2)){
                            ser_level_count[3][Max_level]++;
                          
                            break;
                       }
                       
                       if((f_all[i].rate/f_all[i].sch_times) <=rate_level[3][j]+((rate_level[3][j+1]-rate_level[3][j])/2)){
                            ser_level_count[3][j]++;
                          
                            break;
                       }
                       
                   }
                   for (int j=0;j<=Max_level ;j++){
                       
                       if((f_all[i].utility/f_all[i].sch_times) > utility_level[3][Max_level-1]+((utility_level[3][Max_level]-utility_level[3][Max_level-1])/2)){
                            uti_level_count[3][Max_level]++;
                           
                            break;
                       }
                       
                       if((f_all[i].utility/f_all[i].sch_times) <=utility_level[3][j]+((utility_level[3][j+1]-utility_level[0][j])/2)){
                            uti_level_count[3][j]++;
                            
                            break;
                       }
                       
                   }
               }
           }
           
    } 
    fout2<<endl;
    float total_f=0;
    for(int i=0;i<4;i++){
           ave_rate[i]= ave_rate[i]/ num_count[i];
           ave_utility[i]= ave_utility[i]/ num_count[i];
           fout2<<"Service type: "<<i <<" ave rate: " << ave_rate[i]<<" ave uti: "<<  ave_utility[i] <<endl;
           fout2<<"number of no gch flow: "<< no_gch_f[i]<<endl;
    }
    fout2<<endl;
    fout2<<"ave_utility for all(sum of utility): "<< ave_u/total_f <<"\n"<<endl;
    for(int i=0 ; i<4 ; i++){
            for(int j=0 ; j<=Max_level ;j++){
                    fout2<<"Service: "<< i<<" rate level: "<<rate_level[i][j]<<" number: "<<ser_level_count[i][j]<<endl;
            }
            fout2<<endl;
    }
    fout2<<endl;
    for(int i=0 ; i<4 ; i++){
            for(int j=0 ; j<=Max_level ;j++){
                    fout2<<"Service: "<< i<<" utility level: "<<utility_level[i][j]<<" number: "<<uti_level_count[i][j]<<endl;
            }
            fout2<<endl;
    }
    fout2<<"no channel times: " << SU_no_ch << endl;
    fout2<<"total sch times: "<< s_t <<endl;
    
    
	system("pause");
}





