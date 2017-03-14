#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

struct graph_node {
	unsigned int val;
	struct graph_node* next;
} *row_ptrs_head, *row_ptrs_tail;

double norm2 (double pagerank[], int length) {
	int i = 0;
	double temp = 0;
	for (i = 0; i < length; i++) {
		temp += pagerank[i] * pagerank[i];
	}
	return sqrt(temp);
}

void addarray(double* origin, double* newarray, int len) {
	int i = 0;
	for (i = 0; i < len; i++) {
		if (newarray[i] != 0)
			origin[i] += newarray[i];
	}
}


struct thread_comm {
	int myrank;
	int num_nodes;
	pthread_attr_t attr;
	// int ind;
};

struct timeval starttime, endtime;
double** temp_buffer;
double* pagerank;
int num_threads;
int num_nodes;
int *split_file;
int* col_inds;
double *col_prob;
pthread_barrier_t barrier0, barrier1, barrier2;

int** col_inds_partition;
struct graph_node** row_ptrs_partition;
struct graph_node* tmp_row_ptrs;



void *pagerank_thread (void* arg) {
	int myrank = *(int *)arg;
	int i, j, k;
	int count = 0;
	int partition_part = (num_nodes / num_threads) + 1;
	
	double last_norm = norm2(pagerank, num_nodes);
	
	pthread_barrier_wait(&barrier0);
	
	//////////////////////////////////////////////////////////////
	if (myrank == 0) {
		gettimeofday(&starttime, NULL);
	}
	double temp;
	while(1) {
		/* Read pagerank and compute own part, then store to temp_buffer[i] */
		
		struct graph_node *current = row_ptrs_partition[myrank];
		i = 0;
		while (current->next != NULL) {
			temp = 0;
			for (j = current->val; j < current->next->val; j++) {
				temp += pagerank[col_inds_partition[myrank][j]] * col_prob[col_inds_partition[myrank][j]];
			}
			temp_buffer[myrank][i] = temp;
			current = current->next;
			i++;
		}
		
		pthread_barrier_wait(&barrier1);
		
		/* Aggregate own part of temp_buffer to corresponding part of pagerank */
		
		for (i = (partition_part * myrank); i < (partition_part * (myrank + 1)) && i < num_nodes; i++) {
			pagerank[i] = 0;
		}
		for (i = 0; i < num_threads; i++) {
			for (j = (partition_part * myrank); j < (partition_part * (myrank + 1)) && j < num_nodes; j++) {
				pagerank[j] += temp_buffer[i][j];
			}
		}
		
		pthread_barrier_wait(&barrier2);
		
		/* Check whether to stop */
		count++;
		if (myrank == 0) {
			printf("Current norm: %f\n", last_norm);
		}
		
		if (fabs(norm2(pagerank, num_nodes) - last_norm) < 0.0000001) {
			if (myrank == 0) {
				gettimeofday(&endtime, NULL);
				printf("iterate count: %d\n", count);
			}
			break;
		}
		
		last_norm = norm2(pagerank, num_nodes);
	}
	
	pthread_exit(NULL);
}

int main (int argc, char* argv[]) {
	int i, j, k;
	num_threads = atoi(argv[3]);
	char text[32];
	
	/****************** Read split file **********************/
	num_nodes = 0;
	char* split_filename = argv[2];
	FILE *fp_split = fopen(split_filename, "r");
	while (fscanf(fp_split, "%s", text) != EOF) {
		num_nodes++;
	}
	// printf("num of nodes: %d\n", num_nodes);
	memset(text, 0, 32);
	
	rewind(fp_split);
	
	split_file = (int*)malloc(sizeof(int) * num_nodes);
	i = 0;
	while (fscanf(fp_split, "%s", text) != EOF) {
		split_file[i] = atoi(text);
		i++;
		memset(text, 0, 32);
	}
	
	fclose(fp_split);
	
	/***************** Read graph file *****************/
	int state = 0;
	int col_len = 0;
	col_inds = NULL;
	char* filename = argv[1];
	FILE *fp = fopen(filename, "r");
	
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
	
	pagerank = (double*)malloc(sizeof(double) * num_nodes);
	for (i = 0; i < num_nodes; i++) {
		pagerank[i] = 1.0 / num_nodes;
	}
	
	col_prob = (double*)malloc(sizeof(double) * num_nodes);
	struct graph_node *current = row_ptrs_head;
	i = 0;
	while (current->next != NULL) {
		col_prob[i] = 1 / (double)(current->next->val - current->val);
		current = current->next;
		i++;
	}
	
	/**************************** Partition graph into subgraph *******************************/
	
	col_inds_partition = (int**)malloc(sizeof(int*) * num_threads);
	row_ptrs_partition = (struct graph_node**)malloc(sizeof(struct graph_node*) * num_threads);
	for (i = 0; i < num_threads; i++) {
		col_inds_partition[i] = (int*)malloc(sizeof(int) * col_len);
		row_ptrs_partition[i] = (struct graph_node*)malloc(sizeof(struct graph_node));
		row_ptrs_partition[i]->val = 0;
	}
	
	for (i = 0; i < num_threads; i++) {
		int row_val = 0;
		current = row_ptrs_head;
		tmp_row_ptrs = row_ptrs_partition[i];
		while (current->next != NULL) {
			for (j = current->val; j < current->next->val; j++) {
				if ((((num_nodes / num_threads) + 1) * i) <= col_inds[j] && col_inds[j] < (((num_nodes / num_threads) + 1) * (i + 1))) {
					col_inds_partition[i][row_val] = col_inds[j];
					row_val++;
				}
			}
			
			tmp_row_ptrs->next = (struct graph_node*)malloc(sizeof(struct graph_node));
			tmp_row_ptrs->next->val = row_val;
			tmp_row_ptrs = tmp_row_ptrs->next;
			
			current = current->next;
		}
	}
	
	/************************* Calculate *************************/
	
	if (num_threads == 1) {
		double* temp_pagerank = (double*)malloc(sizeof(double) * num_nodes);
		
		double last_norm = norm2(pagerank, num_nodes);
		printf("original norm: %f\n", last_norm);
		
		int count = 0;
		gettimeofday(&starttime, NULL);
		while (1) {
			current = row_ptrs_head;
			i = 0;
			while (current->next != NULL) {
				double temp = 0;
				for (j = current->val; j < current->next->val; j++) {
					temp += pagerank[col_inds[j]] * col_prob[col_inds[j]];
				}
				temp_pagerank[i] = temp;
				
				current = current->next;
				i++;
			}
			
			count++;
			printf("Current norm: %f\n", norm2(temp_pagerank, num_nodes));
			if (fabs(norm2(temp_pagerank, num_nodes) - last_norm) < 0.0000001) {
				memcpy(pagerank, temp_pagerank, sizeof(double) * num_nodes);
				memset(temp_pagerank, 0, sizeof(double) * num_nodes);
				break;
			}
			
			last_norm = norm2(temp_pagerank, num_nodes);
			memcpy(pagerank, temp_pagerank, sizeof(double) * num_nodes);
			memset(temp_pagerank, 0, sizeof(double) * num_nodes);
		}
		gettimeofday(&endtime, NULL);
		
		printf("iterate count: %d\n", count);
		printf("time used: %lu microseconds\n", (1000000 * (endtime.tv_sec - starttime.tv_sec) + (endtime.tv_usec - starttime.tv_usec)));
	}
	else {
		temp_buffer = (double**)malloc(sizeof(double*) * num_threads);
		for (i = 0; i < num_threads; i++) {
			temp_buffer[i] = (double*)malloc(sizeof(double) * num_nodes);
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
		
		printf("time used: %lu microseconds\n", (1000000 * (endtime.tv_sec - starttime.tv_sec) + (endtime.tv_usec - starttime.tv_usec)));
	}
	
	/************************** Finalize ***************************/
	printf("writing file......\n");
	FILE *fp_write = fopen("./mypagerank.result", "w");
	fprintf(fp_write, "time: %lu ms\n", 1000000 * (endtime.tv_sec - starttime.tv_sec) + (endtime.tv_usec - starttime.tv_usec));
	fprintf(fp_write, "node_id\tpagerank\n");
	for (i = 0; i < num_nodes; i++) {
		fprintf(fp_write, "%d\t\t%.12lf\n", i, pagerank[i]);
	}
	fclose(fp_write);
	printf("Finish writing file\n");
}
