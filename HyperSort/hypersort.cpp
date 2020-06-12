#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <iostream>
#include <limits>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <malloc.h>
#include <string>
#include <string.h>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <bits/stdc++.h>
#include <utility>
#include <cctype>
#include <tuple>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define SIZE 6400000

int nproc;
int count=0;

using namespace std;


void sort(int *a, int size){

	for(int j=1;j<size;j++){
		for(int i=0;i<(size-j);i++){
			if(a[i]>a[i+1]){
				signed int temp=a[i];
				a[i]=a[i+1];
				a[i+1]=temp;
			}
		}
	}

}

void write_output(int num_input, signed int sorted_output[])
{
  std::ofstream outfile;
  outfile.open("outfile5.txt");

  for(int i=0;i<num_input;i++){
  	outfile<<sorted_output[i]<<" ";
  }
  
  outfile.close();
}

double timer()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

int main(int argc, char* argv[])
{
	int curr_proc, num_proc, iproc;
	signed int *unsorted;
	MPI_Status status;
	signed int *process_block; 
	int num_input;
	MPI_Init (&argc, &argv);      /* starts MPI */
  	MPI_Comm_rank (MPI_COMM_WORLD, &curr_proc);        /* get current process id */
  	MPI_Comm_size (MPI_COMM_WORLD, &num_proc);        /* get number of processes */
  	double starttime;

  	string input= "input5.txt" ;

		std::ifstream in(input);
		string line;
		int line_num = 0;

		while (getline(in, line))
		{			
			    if(line_num == 0){
			    	num_input = std::stoi(line);	
			    	unsorted = (int *) malloc (num_input * sizeof (int));
			    //	cout<<num_payload<<endl;	    	
			    }
			    if(line_num == 1){
			    	string str = "";
		    		int array_counter = 0 ;
		    		for(int i=0 ; i < line.length(); i++){  			
		    			char c = line.at(i) ;
		    			if(c == ' '){
		    				signed int n;
		    				if(str.at(0) == '-'){
		    					int len = str.length()-1;
		    					string str1 = str.substr(1, len);
		    					n = std::stoi(str1);
		    					
		    					n = n*(-1);
		    				}
		    				else{
		    					n = std::stoi(str);
		    				}
		    				
		    				unsorted[array_counter] = n;
		    				//cout<<"n= "<<n<<" line "<<array_counter<<endl;
		    				array_counter++;
		    				str=""; 	
		    			}
		    			else{
		    				str += c ;
		    			}
		    		}
		    	//	unsorted[array_counter] = std::stoi(str);
			    }
			line_num++ ;
		}

  	int num_data_proc = num_input/num_proc;
  //	cout<<"num_data_proc "<<num_data_proc;
	process_block = (int *) malloc (num_input * sizeof (int));
	for(int i=0;i<num_data_proc;i++){
		process_block[i]= unsorted[(num_data_proc*curr_proc)+i];
	}

	sort(process_block, num_data_proc);

	MPI_Comm newcomm;
	MPI_Comm_split(MPI_COMM_WORLD, 1, iproc, &newcomm);
	int groupSize= num_proc;
	starttime = timer();
	while(groupSize>1){
		cout<<"#################WHILE LOOP STARTS#####################"<<endl;
		MPI_Comm_size(newcomm, &nproc);
		MPI_Comm_rank(newcomm, &iproc);

		int* pivot=(int*)malloc(sizeof(int));
		*pivot = process_block[num_data_proc/2];
	//	cout<<"PIVOT "<<*pivot<<endl;
		int* recv_pivot=(int*)malloc(sizeof(int));

			if((iproc%groupSize)==0){

				for(int i=1;i<nproc;i++){
					MPI_Send(pivot, 1, MPI_INT, i, 0, newcomm);
				//	cout<<"SENT PIVOT "<<*pivot<<endl;
				}
				*recv_pivot = *pivot;
			}
			else{
			//	MPI_Send(pivot, 1, MPI_INT, 0, 0, newcomm);
				MPI_Recv(recv_pivot, 1, MPI_INT, 0, 0, newcomm, &status);
			//	cout<<"RECEIVED PIVOT "<<*recv_pivot<<endl;
			}
		//}
		
		signed int *sending_buf=(int *) malloc (num_data_proc * sizeof (int));
		int* sending_buf_counter=(int*)malloc(sizeof(int));
		*sending_buf_counter=0;
		signed int *received=(int *) malloc (num_input * sizeof (int));
		int* rec_counter=(int*)malloc(sizeof(int));
		*rec_counter=0;
		for(int i=0; i < num_data_proc;i++){
			if(iproc < nproc/2){
				if(process_block[i] >= *recv_pivot){
					sending_buf[*sending_buf_counter] = process_block[i];
					(*sending_buf_counter)++;
				}
			}
			else{
				if(process_block[i] <= *recv_pivot){
					sending_buf[*sending_buf_counter] = process_block[i];
					(*sending_buf_counter)++;
				}
			}
		}

		if(iproc < nproc/2){
			MPI_Send(sending_buf_counter, 1, MPI_INT, iproc+(nproc/2),0, newcomm);
			MPI_Send(sending_buf, *sending_buf_counter, MPI_INT, iproc+(nproc/2),0, newcomm);
			MPI_Recv(rec_counter, 1, MPI_INT, iproc+(nproc/2),0,newcomm,&status);
			MPI_Recv(received, *rec_counter, MPI_INT, iproc+(nproc/2),0,newcomm,&status);

		}
		else{
			MPI_Recv(rec_counter, 1, MPI_INT, iproc-(nproc/2),0,newcomm,&status);
			MPI_Recv(received, *rec_counter, MPI_INT, iproc-(nproc/2),0,newcomm,&status);
			MPI_Send(sending_buf_counter, 1, MPI_INT, iproc-(nproc/2),0, newcomm);
			MPI_Send(sending_buf, *sending_buf_counter, MPI_INT, iproc-(nproc/2),0, newcomm);

		}

		if(iproc<nproc/2){
			num_data_proc = num_data_proc-*sending_buf_counter;
			for(int i=0;i<*rec_counter;i++){
				process_block[num_data_proc]=received[i];
				num_data_proc++;
			}
		}
		else{
			int offset=0;
			for(int i=0;i<num_data_proc;i++){
				if(process_block[i]>=*recv_pivot){
					received[*rec_counter+offset]=process_block[i];
					offset++;
				}
			}
			signed int *temp = process_block;
			process_block=received;
			received=temp;
			num_data_proc=*rec_counter+(num_data_proc-*sending_buf_counter);
		}


		sort(process_block, num_data_proc);
		MPI_Comm_split(newcomm, iproc < nproc/2 , iproc, &newcomm);
		groupSize /= 2;
	}
		

	int proc, number;
	MPI_Comm_rank(MPI_COMM_WORLD, &proc);
	MPI_Comm_size (MPI_COMM_WORLD, &number); 
	MPI_Barrier(MPI_COMM_WORLD);
	cout<<"---------TIME TAKEN-------------"<<(timer()-starttime)<<endl;
	int* process_block_counter=(int*)malloc(sizeof(int));
	*process_block_counter= num_data_proc;
	signed int *received_proc=(int *) malloc (num_input * sizeof (int));
	int* proc_counter=(int*)malloc(sizeof(int));
	*proc_counter=0;
	std::vector<std::vector<signed int>> blocks ;

	if(proc != 0){
		//cout<<"sending blocks"<<endl;
			MPI_Send(process_block_counter, 1, MPI_INT, 0,0, MPI_COMM_WORLD);
			MPI_Send(process_block, *process_block_counter, MPI_INT, 0, 0, MPI_COMM_WORLD);
		//	cout<<"block for process "<<proc<<" sent"<<endl;
	}else{
		//cout<<"receiving blocks"<<endl;
		for (int i=1; i<number; i++) {
            MPI_Recv(proc_counter, 1, MPI_INT, i,0,MPI_COMM_WORLD,&status);
           // cout<<"proc_counter for process "<<i<<" recieved is "<<*proc_counter<<endl;
			MPI_Recv(received_proc, *proc_counter, MPI_INT, i,0,MPI_COMM_WORLD,&status);
			std::vector<signed int> tmp ;
			//cout<<"From process "<<i<<" : ";
			for (int j = 0; j < *proc_counter; ++j)
				{
				//	cout<<received_proc[j]<<" ";
					tmp.push_back(received_proc[j]);
				}
				//cout<<""<<endl;	
           if(tmp.size() > 0){
				blocks.push_back(tmp);
			}
		}
		std::vector<signed int> tmp1 ;
		for (int j = 0; j < num_data_proc; ++j)
				{
				//	cout<<received_proc[j]<<" ";
					tmp1.push_back(process_block[j]);
				}
				//cout<<""<<endl;
			if(tmp1.size() > 0){
				blocks.push_back(tmp1);
			}	
            


		signed int sorted_output[num_input] ;
		int sorted_output_counter=0;
	//	cout<<"Final list "<<endl;
		while(blocks.size()>0){
			int least = std::numeric_limits<int>::max();
			std::vector<signed int> least_vect;
			int index;
			for(int k=0; k < blocks.size(); k++){
				if(blocks[k][0] < least){
					least = blocks[k][0];
					least_vect = blocks[k];
					index = k;
				}
			}
			//cout<<"least vector"<<endl;
			for(int i=0; i<least_vect.size(); i++){
				sorted_output[sorted_output_counter] = least_vect[i];
				sorted_output_counter++;
			}
			//cout<<""<<endl;
			blocks.erase(blocks.begin() + index);

		}

  		write_output(num_input, sorted_output);
	}
  	
  		MPI_Finalize();
  	return 0;
}
