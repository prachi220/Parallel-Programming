#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> 
#include <pthread.h>
#include <math.h>
#include <omp.h>
#include <openssl/md5.h>

#define SIZE 1024*1024*16

int *arr ;
int *final_arr ;
int d ;
int B;
int *sum ;
int *incr ;
int ARRAY_SIZE ;
int N ;

 /*void *top_down(void * threadId){
 	int i;
 	 i = *((int *) threadId);
  // printf("i is %d\n", (i)) ;
   int power = pow(2,d) ;
  int power1 = pow(2,d+1) ;
 	sum[i + power1 -1] = sum[i + power - 1] + sum[i + power1 -1] ;

 }

 void *down_sweep(void * threadId){
 	int i;
 	i = *((int *) threadId);
  int power = pow(2,d) ;
  int power1 = pow(2,d+1) ;
 	int t = sum[i + power -1] ;
 	sum[i + power -1] = sum[i + power1 - 1] ;
 	sum[i + power1 - 1] = t + sum[i + power1 - 1] ;
 }

 void scan_sum(int index_arr){
//  printf("%s\n", "SCAN SUM STARTS" );
  int r =  log(index_arr)/log(2) ;
 
 	for(d= 0 ; d < r ; d++){
    int power1 = pow(2, (d+1)) ;
   //  printf("pow d is %f\n", power1);
 		int j = 0 ;
 		int index_thread = 0 ;
 		int NUMTHREADS = index_arr/power1 ;
   //  printf("numthreads %d\n", NUMTHREADS);
 		pthread_t threads[NUMTHREADS];
 		while((j * power1) < index_arr){
      int *arg = malloc(sizeof(*arg));
 			int i = j * power1 ;
      *arg = i ;
 			int rc  = pthread_create(&threads[index_thread], NULL, top_down, arg);
			    if (rc){
			//      printf("ERROR; return code from pthread_create() is %d\n", rc);
			      exit(-1);
			    }
			j++ ;
			index_thread++ ;				
 		}

 		for (int k=0; k<NUMTHREADS; k++){
    		pthread_join(threads[k], NULL);
  		}
 	}

 	last = sum[index_arr-1] ;
 	sum[index_arr-1] = 0 ;
 	
 	for(d= r-1 ; d > -1 ; d--){
 		int j = 0 ;
 		int index_thread1 = 0 ;
 		int power1 = pow(2,(d+1)) ;
 		int NUMTHREADS = index_arr/power1 ;
 		pthread_t threads[NUMTHREADS];
 		while((j * power1) < index_arr){
      int *arg = malloc(sizeof(*arg));
 			int i = j * power1 ;
      *arg = i ;
 			int rc  = pthread_create(&threads[index_thread1], NULL, down_sweep, arg);
			    if (rc){
			      printf("ERROR; return code from pthread_create() is %d\n", rc);
			      exit(-1);
			    }
			j++ ;
			index_thread1++ ;				
 		}

 		for (int k=0; k<NUMTHREADS; k++){
    		pthread_join(threads[k], NULL);
  		}
 	}

 	 
 }

 int processor_sum(int index, int* a){
  int sum = 0 ;
  for(int i= index ; i<B+index ; i++){
    sum = sum + a[i] ;
  }
  return sum ;
 }

 void *proc_top_down(void * threadId){
  int i;
   i = *((int *) threadId);
  // printf("i is %d\n", (i)) ;
   int power = pow(2,d) ;
  int power1 = pow(2,d+1) ;
  arr[i + power1 -1] = arr[i + power - 1] + arr[i + power1 -1] ;

 }

 void *proc_down_sweep(void * threadId){
  int i;
  i = *((int *) threadId);
  int power = pow(2,d) ;
  int power1 = pow(2,d+1) ;
  int t = arr[i + power -1] ;
  arr[i + power -1] = arr[i + power1 - 1] ;
  arr[i + power1 - 1] = t + arr[i + power1 - 1] ;
 }

void proc_prefix_sum(int start_index, int prescan, int index_arr){
  int r =  log(B)/log(2) ;
//  printf("start_index is %d\n", start_index);
//  printf("r is %d\n", r);
 
  for(d= 0 ; d < r ; d++){
    int power1 = pow(2, (d+1)) ;
   //  printf("pow d is %f\n", power1);
    int j = 0 ;
    int index_thread = 0 ;
    int NUMTHREADS = B/power1 ;
   //  printf("numthreads %d\n", NUMTHREADS);
    pthread_t threads[NUMTHREADS];
    while((j * power1) < B){
      int *arg = malloc(sizeof(*arg));
      int i = j * power1 + start_index;
      *arg = i ;
      int rc  = pthread_create(&threads[index_thread], NULL, proc_top_down, arg);
          if (rc){
      //      printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
          }
      j++ ;
      index_thread++ ;        
    }

    for (int k=0; k<NUMTHREADS; k++){
        pthread_join(threads[k], NULL);
      }
  }

  arr[start_index+B-1] = prescan ;
  
  for(d= r-1 ; d > -1 ; d--){
    int j = 0 ;
    int index_thread1 = 0 ;
    int power1 = pow(2,(d+1)) ;
    int NUMTHREADS = B/power1 ;
    pthread_t threads[NUMTHREADS];
    while((j * power1) < B){
      int *arg = malloc(sizeof(*arg));
      int i = j * power1 + start_index ;
      *arg = i ;
      int rc  = pthread_create(&threads[index_thread1], NULL, proc_down_sweep, arg);
          if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
          }
      j++ ;
      index_thread1++ ;       
    }

    for (int k=0; k<NUMTHREADS; k++){
        pthread_join(threads[k], NULL);
      }
  }

 }
*/
void *seq_prefix(void * threadId){
  int i;
   i = *((int *) threadId);
   int final = 0 ;
    for (int j=i; j<i+B; j++){
      arr[j] = final ;
      final = final + arr[j] ;

    }
 }

 void *sum_array(void * threadId){
    int i;
    i = *((int *) threadId);
    printf("i %d\n", i);
    int t = i/N ;
    int add = 0 ;
    for(int j=i ; j< i+ N ; j++){
      add = add + arr[j] ;
    }
    sum[t] = add ;
    
 }

 void *final_array(void * threadId){
    int i;
    i = *((int *) threadId);
    int a[N] ;
    int t = i/N ;
    // printf("i %d\n", i);
    // printf("t %d\n", t);
    // printf("incr %d\n", incr[t]);
    a[0] = incr[t] ;
    for(int j=1 ; j<N; j++){
      a[j] = arr[j+i-1] ;
    }

    // for(int j=0 ; j< N; j++){
    //   printf("a = %d\n", a[j]);
    // }

    final_arr[i] = a[0] ;
    for (int j=1; j<N; j++){
      final_arr[j+i] = final_arr[j+i-1] + a[j] ;
    }

    // for(int j=i ; j< N+i; j++){
    //   printf("j %d", j);
    //   printf("final_arr = %d\n", final_arr[j]);
    // } 
 }
 

int main(int argc, char *argv[]) {

  int threads = 0 ;

  if(argc == 2){
    threads = atoi(argv[1]) ;
  }
  else{
    printf("%s\n", "Too many arguments passed");
  }

  double  dif;
 
  arr = malloc (sizeof(int)*SIZE);


  FILE* file = fopen("input.txt", "r"); /* should check the result */
  char c;
  char buffer[256] ;
  int index_buf = 0;
 // printf("%s\n","malloc successful" );

  while ((c = getc(file)) != EOF)
    {
        
   // printf("%s",  "hey there");
      if(c == '\n'){
        ARRAY_SIZE = atoi(buffer) ;
        memset(buffer, 0, sizeof(buffer));
        index_buf = 0 ;
        break ;
      }
      else{
        buffer[index_buf] = c ;
        index_buf++ ;
   //     printf("%c\n", buffer[index_buf]);
        
      }
      
    }
    printf("Array: %d\n", ARRAY_SIZE);

    int *temp = realloc(arr, ARRAY_SIZE*sizeof(int));
    if ( temp != NULL ) //realloc was successful
    {
    arr = temp;
    }
    else //there was an error
    {
    free(arr);
    printf("Error allocating memory!\n");
    return 1;
    }

   //  printf("%s\n","realloc successful" );
   
     int index_arr1 = 0 ;

      while ((c = getc(file)) != EOF)
    {
        
   // printf("%s",  "hey there");
    //  printf("c= %c\n", c );

        if(c == ' '){
          int num = atoi(buffer) ;
        //  printf("num : %d\n", num);
          arr[index_arr1] = num ; 
          index_arr1++ ;
          memset(buffer, 0, sizeof(buffer));
          index_buf = 0 ;
        }
        else{
          buffer[index_buf] = c ;
          index_buf++ ;
      }
      
    }
     int num1 = atoi(buffer) ;
     printf("num : %d\n", num1);
     arr[index_arr1] = num1 ; 
     index_arr1++ ;
    
    // printf("index : %d\n", index_buf);

    // printf("%d\n", ARRAY_SIZE);
    // for(int y=0 ; y< ARRAY_SIZE ; y++){
    //   printf("y: %d\n", y);
    //   printf("arr : %d", arr[y] );
    // }

  fclose(file);
   B= threads ;
   N = ARRAY_SIZE/B ;

   final_arr = malloc(sizeof(int)*(ARRAY_SIZE+1)) ;
  sum = malloc (sizeof(int)*B);
  incr = malloc (sizeof(int)*B);

//CODE HERE
    double start = omp_get_wtime();

    int index_thread1 = 0 ;

    pthread_t Threads[B];
    int i = 0 ;
   
    while(i< ARRAY_SIZE){
  //    printf("i = %d\n", i);
      int *arg = malloc(sizeof(*arg));
      *arg = i;
        int rc  = pthread_create(&Threads[index_thread1], NULL, sum_array, arg);
          if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
          }
          index_thread1++ ; 
          i = i + N ;      
    }

    for (int k=0; k<B; k++){
        pthread_join(Threads[k], NULL);
      }

    // for(int i=0 ; i< B ; i++){
    //   printf("sum %d\n", sum[i] );
    // }
    int last = 0 ;
    for(int p=0 ; p < B ; p++){
      last = last + sum[p] ;
    }
     
     incr[0] = 0 ;
    for (int j=1; j<B; j++){
      incr[j] = incr[j-1] + sum[j-1] ;
    }

    // for(int i=0 ; i< B ; i++){
    //   printf("prscan %d\n", incr[i] );
    // }

    int i1 = 0 ;
    int index_thread = 0 ;
    pthread_t Threads1[B];
    while(i1< ARRAY_SIZE){
  //    printf("i = %d\n", i);
      int *arg = malloc(sizeof(*arg));
      *arg = i1;
        int rc  = pthread_create(&Threads1[index_thread], NULL, final_array, arg);
          if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
          }
          index_thread++ ; 
          i1 = i1 + N ;      
    }

    for (int k=0; k<B; k++){
        pthread_join(Threads1[k], NULL);
      }
      final_arr[ARRAY_SIZE] = last ;
    for(int j=0 ; j< ARRAY_SIZE+1; j++){
      printf("final_arr = %d\n", final_arr[j]);
    } 

    
 //  // prefix_sum(arr, index_arr) ;
  double end = omp_get_wtime();// end the timer
  dif = end - start ;
  printf("difference %f\n", dif );

 

  FILE *f1 = fopen("temp1.txt", "w");
  if (f1 == NULL)
  {
      printf("Error opening file!\n");
      exit(1);
  }
  for(int i=1 ; i< ARRAY_SIZE+1 ; i++){
    fprintf(f1, "%d", final_arr[i]);
    fprintf(f1, "%s", " ");
  }

  fclose(f1) ;

  char *line ;
  FILE* file2 = fopen("temp1.txt", "r");

   while (fgets(line, sizeof(line), file2)) {
     
  }

  fclose(file2) ;

   FILE *f = fopen("output.txt", "w");
  if (f == NULL)
  {
      printf("Error opening file!\n");
      exit(1);
  }

  unsigned char digest[16];
 // const char* string = "Hello World";
  MD5_CTX context;
  MD5_Init(&context);
  MD5_Update(&context, line, strlen(line));
  MD5_Final(digest, &context);

  char *out = (char*)malloc(33);
  
  for (int n = 0; n < 16; ++n) {
      snprintf(&(out[n*2]), 16*2, "%02x", (unsigned int)digest[n]);
  }

  printf("%s\n", out );


  /* print some text */
  const char *text = "Threads: ";
  fprintf(f, "%s", text);
  fprintf(f, "%d\n", threads) ;

  const char *text1 = "Time: ";
  fprintf(f, "%s", text1);
  fprintf(f, "%f\n", dif) ;

  const char *text2 = "Md5-sum: ";
  fprintf(f, "%s", text2);
  fprintf(f, "%s\n", out );

  printf("%s\n", "hello");

  fclose(f);
  
    
}