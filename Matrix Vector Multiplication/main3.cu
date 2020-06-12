#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>


#define WARP_SIZE 32
#define BLOCK_SIZE 1024


#define gpuErrChk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess)
    {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
    }
}


__global__  void matrixMultiply(int num_rows, int * ptr , int * indices, unsigned long int * data , unsigned long int * x, unsigned long int * y){
	 

//	int  row = blockDim.x * blockIdx.x + threadIdx.x;
//	if(row < num_rows ){
//		unsigned long int dot = 0;
//		int  row_start = ptr[row];
//		int  row_end    = ptr[row +1];
//		for (int jj = row_start; jj < row_end; jj++)
//			dot += data[jj] * x[indices[jj]];
//		y[row] += dot;
//	}

	 	__shared__ unsigned long int sdata[1024];
     __shared__ int ptrs[1024/WARP_SIZE][2];
    
     const int thread_id   = BLOCK_SIZE * blockIdx.x + threadIdx.x;  // global thread index
     const int thread_lane = threadIdx.x & (WARP_SIZE-1);            // thread index within the warp
     const int warp_id     = thread_id   / WARP_SIZE;                // global warp index
     const int warp_lane   = threadIdx.x / WARP_SIZE;                // warp index within the CTA
     const int num_warps   = (BLOCK_SIZE / WARP_SIZE) * gridDim.x;   // total number of active warps

     for(int row = warp_id; row < num_rows; row += num_warps){
         // use two threads to fetch ptr[row] and ptr[row+1]
         // this is considerably faster than the more straightforward option
         if(thread_lane < 2)
             ptrs[warp_lane][thread_lane] = ptr[row + thread_lane];
         const int row_start = ptrs[warp_lane][0]; //same as: row_start = ptr[row];
         const int row_end   = ptrs[warp_lane][1]; //same as: row_end   = ptr[row+1];

         // compute local sum
         sdata[threadIdx.x] = 0;
         for(int jj = row_start + thread_lane; jj < row_end; jj += WARP_SIZE)
             sdata[threadIdx.x] += data[jj] * x[indices[jj]];

         // reduce local sums to row sum (ASSUME: warpsize 32)
         if (thread_lane < 16) { sdata[threadIdx.x] += sdata[threadIdx.x + 16]; __syncthreads(); }
         if (thread_lane <  8) { sdata[threadIdx.x] += sdata[threadIdx.x +  8]; __syncthreads(); }
         if (thread_lane <  4) { sdata[threadIdx.x] += sdata[threadIdx.x +  4]; __syncthreads(); }
         if (thread_lane <  2) { sdata[threadIdx.x] += sdata[threadIdx.x +  2]; __syncthreads(); }
         if (thread_lane <  1) { sdata[threadIdx.x] += sdata[threadIdx.x +  1]; __syncthreads(); }

         // first thread writes warp result
         if (thread_lane == 0)
             y[row] += sdata[threadIdx.x];
     }
}


int main (int argc, char **argv)
{
    int rank, nprocs;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);

    FILE* file = fopen(argv[1], "r"); /* should check the result */
  	char c;
  	char buffer[256] ;
  	int index_buf = 0;
  	int line_num=0;
  	int take_input=0;
  	int take_input_A=0;
  	int take_input_B=0;
  	int ARRAY_SIZE;
  	int num_rows;
  	unsigned long int *B, *device_B;
  	int cpu_counter=0, num_non_zero=0;
  	unsigned long int *result_host, *result_gpu;
  	//int result;
  	int Acol=0, Arow=0;
  	int num, line_num_B=0;
  	int  totalBlocks;
  	unsigned long int *data, *device_data;
  	int *indices, *rows, *ptr, *device_indices, *device_rows, *device_ptr;

  	while ((c = getc(file)) != EOF)
    {
    	//printf("c= %c\n", c );
    	if(line_num==0){
    		if(c=='\n'){
				line_num++;
			}
    	}
    	if(line_num == 1){
		   	if(c == ' '){
		   		take_input= 1;
		   	}
		   	if(take_input){
		   		if(c == '\n'){
		   			ARRAY_SIZE = atoi(buffer) ;
		   			
		   			num_rows=ARRAY_SIZE/nprocs;
				if(rank==(nprocs-1)){
					num_rows = num_rows+(ARRAY_SIZE%nprocs);
				}
        			memset(buffer, 0, sizeof(buffer));
        			index_buf = 0 ;
        			line_num++;
        			B = (unsigned long int  *)malloc (sizeof(unsigned long int )*ARRAY_SIZE);
        			data = (unsigned long int  *)malloc (sizeof(unsigned long int )*ARRAY_SIZE*num_rows);
        			indices = (int *)malloc (sizeof(int)*ARRAY_SIZE*num_rows);
        			rows = (int *)malloc (sizeof(int)*ARRAY_SIZE*num_rows);
        			ptr = (int *)malloc (sizeof(int)*(num_rows+1));

		   		}
		   		else{
		   			buffer[index_buf] = c ;
		   			index_buf++;
		   		}
		   	}	    	
		}
		else{

			if(line_num > 1){

				if(c == 'A'){
					take_input_A=1;
					take_input=0;
				}
				if(take_input_A){
					if(line_num==2 && (c=='\n')){
						line_num++;
					}else{
					//printf("line_num: %d\n", line_num);
						if(c == ' '){
							
							num = atoi(buffer) ;
							
							if(Acol==0){
								Arow=num;
							}
							// printf("rank: %d", rank);
							// printf(" Arow: %d", Arow);
							// printf(" Acol: %d", Acol);
							// printf(" num: %d\n", num);
							if((Arow >= (ARRAY_SIZE/nprocs)*rank) && (Arow < (((ARRAY_SIZE/nprocs)*rank)+num_rows))){
								if(Acol==0){
									rows[num_non_zero]=Arow-(rank*(ARRAY_SIZE/nprocs));
								}
								if(Acol==1){
									indices[num_non_zero]=num;
								}
							}
							Acol++;
							memset(buffer, 0, sizeof(buffer));
	    					index_buf = 0 ;
						}else{
							if(c== '\n'){
								unsigned long int num1 = strtoul(buffer, NULL, 0) ;
								// printf("rank: %d", rank);
								// printf(" Arow: %d", Arow);
								// printf(" Acol: %d", Acol);
								// printf(" num: %d\n", num);
								if((Arow >= (ARRAY_SIZE/nprocs)*rank) && (Arow < (((ARRAY_SIZE/nprocs)*rank)+num_rows))){
									data[num_non_zero]=num1;
									num_non_zero++;
								}
								Acol=0;
								memset(buffer, 0, sizeof(buffer));
		    					index_buf = 0 ;
		    					line_num++;
							}else{
								if(c != 'A'){
									buffer[index_buf] = c;
									index_buf++;
								}
							}
						}
					}
				}

				if(c == 'B'){
					take_input_B=1;
					take_input_A=0;
					line_num_B=line_num;
					index_buf=0;
				}
				if(take_input_B){
					// printf("c= %c\n", c );
					if(c== '\n'){
						if(line_num==line_num_B){
							line_num++;
						}else{
							unsigned long int num1 = strtoul(buffer, NULL, 0) ;
							B[line_num-line_num_B-1]=num1;
	        				memset(buffer, 0, sizeof(buffer));
	        				index_buf = 0 ;
	       					line_num++;
						}		
					}
					else{
						if(c != 'B'){
								buffer[index_buf] = c;
								index_buf++;
						}	
					}
				}
			}
		}
	}
	unsigned long int num1 = strtoul(buffer, NULL, 0) ;
	//printf(" num: %d\n", num);
	B[line_num-line_num_B-1]=num1;
	memset(buffer, 0, sizeof(buffer));
	index_buf = 0 ;
	fclose(file);

	//printf(" Printing data here: ");
	//printf("rank: %d", rank);
	//for(int i=0;i<num_non_zero;i++){
	//	printf(" %d", data[i] );
	//}
	//printf("\n");

	int count=0;
	int ptr_counter=1;
	ptr[0]=0;
	for(int n=0; n < num_rows; n++){
		for(int i=0;i<num_non_zero;i++){
			if(n==rows[i]){
				count++;
			}
		}
		ptr[ptr_counter]= ptr[ptr_counter-1]+count;
		ptr_counter++;
		count=0;
	}


	result_host = (unsigned long int *)malloc(sizeof(unsigned long int) * num_rows);
    //GPU

	gpuErrChk(cudaMalloc((void**)&device_data, sizeof(unsigned long int ) * ARRAY_SIZE * num_rows));
    gpuErrChk(cudaMalloc((void**)&device_B, sizeof(unsigned long int ) * ARRAY_SIZE));
    gpuErrChk(cudaMalloc((void**)&device_indices, sizeof(int) * ARRAY_SIZE * num_rows));
    gpuErrChk(cudaMalloc((void**)&device_ptr, sizeof(int) * (num_rows+1)));
	gpuErrChk(cudaMalloc((void**)&result_gpu, sizeof(unsigned long int) * num_rows));
	
	gpuErrChk(cudaMemcpy(device_data, data, sizeof(unsigned long int )*ARRAY_SIZE*num_rows, cudaMemcpyHostToDevice));
	gpuErrChk(cudaMemcpy(device_B, B, sizeof(unsigned long int ) * ARRAY_SIZE, cudaMemcpyHostToDevice));
	gpuErrChk(cudaMemcpy(device_indices, indices, sizeof(int) * ARRAY_SIZE * num_rows, cudaMemcpyHostToDevice));
	gpuErrChk(cudaMemcpy(device_ptr, ptr, sizeof(int) * (num_rows+1), cudaMemcpyHostToDevice));
	//10,1024
	//int non_zero_rows=ptr_counter-1;

	memset(result_host, 0, sizeof(result_host));
	cudaMemset(result_gpu, 0, sizeof(unsigned long int) * num_rows);
	int NUM_BLOCKS = ((ARRAY_SIZE * num_rows)/1024) + 1;
	matrixMultiply<<<NUM_BLOCKS, 1024>>>(num_rows, device_ptr, device_indices,device_data, device_B, result_gpu);

    cudaThreadSynchronize(); 
    

    //@@ Copy the GPU memory back to the CPU here
    gpuErrChk(cudaMemcpy(result_host, result_gpu, sizeof(unsigned long int) * num_rows, cudaMemcpyDeviceToHost));
     

    //@@ Free the GPU memory here
    //cudaFree(device_data);
    cudaFree(device_B);
    cudaFree(device_indices);
    cudaFree(device_ptr);
    cudaFree(result_gpu);

   // printf(" gather_count : %d\n",*gather_count);
	//  printf(" num rows : %d\n",num_rows);
	//for(int i=0;i<num_rows;i++){
	//	printf("rank : %d value of result : %ld \n", rank, result_host[i] );
	//}
	//printf("\n"); 

	MPI_Barrier(MPI_COMM_WORLD);
	unsigned long int * final_answer = (unsigned long int *)malloc(sizeof(unsigned long int) * ARRAY_SIZE);
	int *gather_count = (int*)malloc(sizeof(int));
	
	//num_rows=ARRAY_SIZE/nprocs;
	//if(rank==(nprocs-1)){
	//	num_rows = num_rows+ARRAY_SIZE%nprocs;
	//}
	*gather_count=ARRAY_SIZE/nprocs; 
	int last= nprocs-1;

	
	
	MPI_Gather(result_host, *gather_count, MPI_UNSIGNED_LONG, final_answer, *gather_count, MPI_UNSIGNED_LONG, last, MPI_COMM_WORLD);

	if(rank==last){
		int remaining = ARRAY_SIZE%nprocs;
		// printf(" remaining : %d\n",remaining);
		for(int i= remaining; i> 0; i--){
		//printf("rank : %llu value of answer : %ld \n", rank, final_answer[i] );
			final_answer[ARRAY_SIZE-i]=result_host[num_rows-i];
		}
	}

	if(rank==last){
		//printf(" Printing final result here");
		//for(int i=0;i<ARRAY_SIZE;i++){
		//	printf("rank : %llu value of answer : %ld \n", rank, final_answer[i] );
		//}

		FILE *f1 = fopen(argv[2], "wb");
	  	if (f1 == NULL)
	  	{
	      printf("Error opening file!\n");
	      exit(1);
	  	}
	  	for(int i=0 ; i< ARRAY_SIZE ; i++){
	    	fprintf(f1, "%ld\n", final_answer[i]);
	  	}

	  	fclose(f1) ;
  	}

  	free(B);
	free(data);
	free(indices);
	free(rows);
	free(ptr);
	free(result_host);
	free(final_answer);

    MPI_Finalize();
    return 0;
    
}

