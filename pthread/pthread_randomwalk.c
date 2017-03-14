#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

struct graph_node {
	unsigned int val;
	struct graph_node* next;
} *row_ptrs_head, *row_ptrs_tail;

int cmp(const void* a, const void *b, void* array) {
	int ia, ib;
	ia = *(int *)a;
	ib = *(int *)b;
	if (((int *)array)[ia] > ((int *)array)[ib])
		return -1;
	else if (((int *)array)[ia] < ((int *)array)[ib])
		return 1;
	else
		return 0;
}

int* getrank(int* input, int len) {
	int i, j;
	int* temp = (int*)malloc(sizeof(int) * len);
	for (i = 0; i < len; i++) {
		temp[i] = i;
	}
	
	qsort_r(temp, len, sizeof(int), cmp, (void *)input);
	
	return temp;
}

struct timeval starttime, endtime;
int** graph;
int* row_count;
int** temp_buffer;
int* pagerank;
int num_threads;
int num_iteration;
int num_nodes;
pthread_barrier_t barrier0, barrier1, barrier2;

void *pagerank_thread (void* arg) {
	unsigned int seed = time(NULL);
	int myrank = *(int *)arg;
	int i, j, k;
	int temp;
	int partition_part = (num_nodes / num_threads) + 1;
	int* current_buffer = (int *)malloc(sizeof(int) * partition_part);
	memset(current_buffer, 0, sizeof(int) * partition_part);
	
	temp = 0;
	for (j = myrank * partition_part; j < ((myrank + 1) * partition_part) && j < num_nodes; j++) {
		while (pagerank[j] > 0) {
			current_buffer[temp] += 1;
			pagerank[j] -= 1;
		}
		temp += 1;
	}
	
	pthread_barrier_wait(&barrier0);
	
	////////////////////////////////////////////////////////////
	if (myrank == 0) {
		gettimeofday(&starttime, NULL);
	}
	
	for (i = 0; i < num_iteration; i++) {
		for (j =0; j < temp; j++) {
			current_buffer[j] = graph[current_buffer[j]][rand_r(&seed) % row_count[current_buffer[j]]];
		}
		if (myrank == 0) {
			printf("\b\b\b\b\b\b\b\b\b\b\b\b");
			printf("iterate: %d", i);
			fflush(stdout);
		}
	}
	
	for (j = 0; j < temp; j++) {
		temp_buffer[myrank][current_buffer[j]] += 1;
	}
	
	pthread_barrier_wait(&barrier1);
	
	for (i = 0; i < num_threads; i++) {
		for (j = (partition_part * myrank); j < (partition_part * (myrank + 1)) && j < num_nodes; j++) {
			pagerank[j] += temp_buffer[i][j];
		}
	}
	
	if (myrank == 0) {
		gettimeofday(&endtime, NULL);
	}
	
	pthread_exit(NULL);
}

int main (int argc, char* argv[]) {
	if (argc != 6) {
		printf("\n    WRONG ARGUMENTS!!! Please use following format:\n");
		printf("\tpthread_randomwalk [graph] [graph_split] [num_threads] [k] [s]\n\n");
		exit(-1);
	}
	
	num_threads = atoi(argv[3]);
	num_iteration = atoi(argv[4]);
	if (num_threads < 1) {
		printf("Too few threads!\n");
		exit(-1);
	}
	if (num_iteration < 2) {
		printf("Too few iteration!\n");
		exit(-1);
	}
	
	int i, j, k;
	int num_output = atoi(argv[5]);
	char text[32];
	
	/******************** Read split file ************************/
	num_nodes = 0;
	char* split_filename = argv[2];
	FILE *fp_split = fopen(split_filename, "r");
	if (fp_split == NULL) {
		printf("Can't open split file!\n");
		exit(-1);
	}
	
	while (fscanf(fp_split, "%s", text) != EOF) {
		num_nodes++;
	}
	memset(text, 0, 32);
	
	fclose(fp_split);
	
	/********************* Read graph file **********************/
	int state = 0;
	int col_len = 0;
	int* col_inds = NULL;
	char* filename = argv[1];
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("Can't open graph file!\n");
		exit(-1);
	}
	
	memset(text, 0, 32);
	while(fscanf(fp, "%s", text) != EOF) {
		if (strcmp(text, "vals:") == 0) {
			state = 1;
		}
		else if (strcmp(text, "col_inds:") == 0) {
			state = 2;
			col_inds = (int*)malloc(sizeof(int) * col_len);
			i = 0;
		}
		else if (strcmp(text, "row_ptrs:") == 0) {
			state = 3;
			
		}
		else {
			if (state == 1) {
				col_len++;
			}
			else if (state == 2) {
				col_inds[i] = atoi(text);
				i++;
			}
			else if (state == 3) {
				if (row_ptrs_head == NULL) {
					row_ptrs_tail = (struct graph_node*)malloc(sizeof(struct graph_node));
					row_ptrs_tail->val = atoi(text);
					row_ptrs_tail->next = NULL;
					row_ptrs_head = row_ptrs_tail;
				}
				else {
					row_ptrs_tail->next = (struct graph_node*)malloc(sizeof(struct graph_node));
					row_ptrs_tail->next->val = atoi(text);
					row_ptrs_tail->next->next = NULL;
					row_ptrs_tail = row_ptrs_tail->next;
				}
			}
			else
				exit(1);
		}
		
		memset(text, 0, 32);
    }
	
	fclose(fp);
	
	/***************************** Initialize *************************************/
	
	pagerank = (int*)malloc(sizeof(int) * num_nodes);
	for (i = 0; i < num_nodes; i++) {
		pagerank[i] = 1;
	}
	
	graph = (int**)malloc(sizeof(int*) * num_nodes);
	row_count = (int*)malloc(sizeof(int) * num_nodes);
	
	struct graph_node *current = row_ptrs_head, *discard;
	discard = current;
	i = 0;
	while (current->next != NULL) {
		row_count[i] = current->next->val - current->val;
		graph[i] = (int *)malloc(sizeof(int) * row_count[i]);
		for (j = current->val; j < current->next->val; j++) {
			graph[i][j - current->val] = col_inds[j];
		}
		
		current = current->next;
		free(discard);
		discard = current;
		i++;
	}
	
	/**************************** Calculate ******************************/
	
	if (num_threads == 1) {
		time_t t;
		srand((unsigned) time(&t));
		
		int* temp_pagerank = (int*)malloc(sizeof(int) * num_nodes);
		memset(temp_pagerank, 0, sizeof(int) * num_nodes);
		
		gettimeofday(&starttime, NULL);
		
		for (i = 0; i < num_iteration; i++) {
			for (j = 0; j < num_nodes; j++) {
				while (pagerank[j] > 0) {
					temp_pagerank[graph[j][rand() % row_count[j]]] += 1;
					pagerank[j] -= 1;
				}
			}
			
			memcpy(pagerank, temp_pagerank, sizeof(int) * num_nodes);
			memset(temp_pagerank, 0, sizeof(int) * num_nodes);
			printf("\b\b\b\b\b\b\b\b\b\b\b");
			printf("iterate: %d", i);
			fflush(stdout);
		}
		
		gettimeofday(&endtime, NULL);
	}
	else {
		temp_buffer = (int**)malloc(sizeof(int*) * num_threads);
		for (i = 0; i < num_threads; i++) {
			temp_buffer[i] = (int*)malloc(sizeof(int) * num_nodes);
			memset(temp_buffer[i], 0, sizeof(int) * num_nodes);
		}
		
		pthread_barrier_init(&barrier0, NULL, num_threads);
		pthread_barrier_init(&barrier1, NULL, num_threads);
		pthread_barrier_init(&barrier2, NULL, num_threads);
		
		pthread_t thread_id[num_threads];
		int rank[num_threads];
		
		for (i = 0; i < num_threads; i++) {
			rank[i] = i;
			if (pthread_create(&thread_id[i], NULL, pagerank_thread, (void *)&rank[i]) != 0) {
				perror("Can't create new thread.\n");
				exit(-1);
			}
		}
		
		for(i = 0; i < num_threads; i++) {
			pthread_join(thread_id[i], NULL);
		}
	}
	
	printf("\b\b\b\b\b\b\b\b\b\b\b           \nNumber of iteration: %d\n", num_iteration);
	printf("Time used: %lu microseconds\n", (1000000 * (endtime.tv_sec - starttime.tv_sec) + (endtime.tv_usec - starttime.tv_usec)));
		
	int count = 0;
	for (i = 0; i < num_nodes; i++) {
		count += pagerank[i];
	}
	printf("Number of walkers: before %d, after %d\n", num_nodes, count);
	/************************** Finalize ***************************/
	
	int* sorted_id = getrank(pagerank, num_nodes);
	printf("\nResult: Top %d visited nodes (descending order)\n-------------------\n  node\t|  visits\n--------|----------\n", num_output);
	for (i = 0; i < num_output; i++) {
		printf("  %d\t|  %d\n", sorted_id[i], pagerank[sorted_id[i]]);
	}
	printf("-------------------\n");
	
	printf("\nWriting file......\n");
	FILE *fp_write = fopen("./randomwalk.result", "w");
	fprintf(fp_write, "Number of iteration: %d\n", num_iteration);
	fprintf(fp_write, "time: %lu ms\n", 1000000 * (endtime.tv_sec - starttime.tv_sec) + (endtime.tv_usec - starttime.tv_usec));
	fprintf(fp_write, "node_id (top %d visited nodes in descending order)\n", num_output);
	for (i = 0; i < num_output; i++) {
		fprintf(fp_write, "%d\n", sorted_id[i]);
	}
	fclose(fp_write);
	
	printf("\nFinish writing file!\n\n");
}
