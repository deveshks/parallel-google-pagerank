#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

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

double *indice;

int cmp(const void* a, const void *b) {
	int ia, ib;
	ia = *(int *)a;
	ib = *(int *)b;
	if (indice[ia] > indice[ib])
		return 1;
	else if (indice[ia] < indice[ib])
		return -1;
	else
		return 0;
}

int* getrank(double* input, int len) {
	int i, j;
	int* rank = (int*)malloc(sizeof(int) * len);
	int* temp = (int*)malloc(sizeof(int) * len);
	for (i = 0; i < len; i++) {
		temp[i] = i;
	}
	indice = input;
	qsort(temp, len, sizeof(int), cmp);
	
	for (i = 0; i < len; i++) {
		rank[temp[i]] = i;
	}
	return rank;
}

int main (int argc, char* argv[]) {
	struct timeval starttime, endtime;
	int numprocs, myrank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	int i, j, k;
	char text[32];
	
	/****************** Read split file **********************/
	int num_nodes = 0;
	char* split_filename = argv[2];
	FILE *fp_split = fopen(split_filename, "r");
	while (fscanf(fp_split, "%s", text) != EOF) {
		num_nodes++;
	}
	memset(text, 0, 32);
	
	rewind(fp_split);
	
	int *split_file = (int*)malloc(sizeof(int) * num_nodes);
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
	int* col_inds = NULL;
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
	
	/**************************** Initialize **************************/

	struct graph_node *current = row_ptrs_head;
	double *col_prob = (double*)malloc(sizeof(double) * num_nodes);
	i = 0;
	while (current->next != NULL) {
		col_prob[i] = 1 / (double)(current->next->val - current->val);
		current = current->next;
		i++;
	}
	
	/************************* Calculate *************************/
	int tag = 1;
	double last_norm;
	double* pagerank = (double*)malloc(sizeof(double) * num_nodes);
	memset(pagerank, 0, sizeof(double) * num_nodes);
	double* temp_pagerank = (double*)malloc(sizeof(double) * num_nodes);
	memset(temp_pagerank, 0, sizeof(double) * num_nodes);
	double* recv_buf = (double*)malloc(sizeof(double) * num_nodes);
	memset((char*)recv_buf, 0, sizeof(double) * num_nodes);
	
	/******************* Set up mapping rule (compression in transmission) ********************/
	int* pre_count = (int*)malloc(sizeof(int) * numprocs);
	int** pre_mapping = (int**)malloc(sizeof(int*) * numprocs);
	for (i = 0; i < numprocs; i++) {
		pre_mapping[i] = (int *)malloc(sizeof(int) * num_nodes);
	}
	int* post_count = (int*)malloc(sizeof(int) * numprocs);
	int** post_mapping = (int**)malloc(sizeof(int*) * numprocs);
	for (i = 0; i < numprocs; i++) {
		post_mapping[i] = (int *)malloc(sizeof(int) * num_nodes);
	}
	
	for (i = 0; i < numprocs; i++) {
		k = 0;
		for (j = 0; j < num_nodes; j++) {
			if (split_file[j] == i) {
				pre_mapping[i][k] = j;
				k++;
			}
		}
		pre_count[i] = k;
	}
	
	for (i = 0; i < numprocs; i++) {
		int m = 0;
		int c = 0;
		int n = 0;
		current = row_ptrs_head;
		while (current->next != NULL) {
			for (j = current->val; j < current->next->val; j++) {
				if (split_file[col_inds[j]] == i) {
					post_mapping[i][c] = m;
					c++;
					break;
				}
			}
			
			current = current->next;
			m++;
		}
		post_count[i] = c;
	}
	
	/*********************** partition graph into subgraph **********************/
	int* col_inds_partition = (int*)malloc(sizeof(int) * col_len);
	struct graph_node* row_ptrs_partition;
	struct graph_node* tmp_row_ptrs;
	int row_val = 0;
	
	current = row_ptrs_head;
	tmp_row_ptrs = (struct graph_node*)malloc(sizeof(struct graph_node));
	tmp_row_ptrs->val = row_val;
	row_ptrs_partition = tmp_row_ptrs;
	while (current->next != NULL) {
		for (j = current->val; j < current->next->val; j++) {
			if (split_file[col_inds[j]] == myrank) {
				col_inds_partition[row_val] = col_inds[j];
				row_val++;
			}
		}
		
		tmp_row_ptrs->next = (struct graph_node*)malloc(sizeof(struct graph_node));
		tmp_row_ptrs->next->val = row_val;
		tmp_row_ptrs = tmp_row_ptrs->next;
		
		current = current->next;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	/** Master **/
	if (myrank == 0) {		
		int data_recved;
		MPI_Status status;
		
		/** initialize pagerank **/
		for (i = 0; i < num_nodes; i++) {
			pagerank[i] = 1.0 / num_nodes;
		}
		last_norm = norm2(pagerank, num_nodes);
		printf("%d: original norm: %f\n", myrank, last_norm);
		
		int count = 0;
		/** send and recv **/
		gettimeofday(&starttime, NULL);
		while(1) {
			for (i = 1; i < numprocs; i++) {
				for (j = 0; j < pre_count[i]; j++) {
					temp_pagerank[j] = pagerank[pre_mapping[i][j]];
				}
				MPI_Send(temp_pagerank, pre_count[i], MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
				// memset((char*)temp_pagerank, 0, sizeof(double) * num_nodes);
			}
			
			
			
			current = row_ptrs_partition;
			i = 0;
			while (current->next != NULL) {
				double temp = 0;
				for (j = current->val; j < current->next->val; j++) {
					temp += pagerank[col_inds_partition[j]] * col_prob[col_inds_partition[j]];
				}
				temp_pagerank[i] = temp; // original length
				i++;
				
				current = current->next;
			}
			
			
			
			for (i = 1; i < numprocs; i++) {
				MPI_Recv(recv_buf, num_nodes, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
				for (j = 0; j < post_count[status.MPI_SOURCE]; j++) {
					temp_pagerank[post_mapping[status.MPI_SOURCE][j]] += recv_buf[j];
				}
				memset((char*)recv_buf, 0, sizeof(double) * num_nodes);
			}
			
			count++;
			printf("Current norm: %f\n", norm2(temp_pagerank, num_nodes));
			if (fabs(norm2(temp_pagerank, num_nodes) - last_norm) < 0.00001) {
				gettimeofday(&endtime, NULL);
				for (i = 1; i < numprocs; i++) {
					MPI_Send(pagerank, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
				}
				memcpy(pagerank, temp_pagerank, sizeof(double) * num_nodes);
				// memset(temp_pagerank, 0, sizeof(double) * num_nodes);
				break;
			}
			
			last_norm = norm2(temp_pagerank, num_nodes);
			memcpy(pagerank, temp_pagerank, sizeof(double) * num_nodes);
			// memset(temp_pagerank, 0, sizeof(double) * num_nodes);
		}
		
		printf("iterate count: %d\n", count);
		printf("time used: %lu microseconds\n", (1000000 * (endtime.tv_sec - starttime.tv_sec) + (endtime.tv_usec - starttime.tv_usec)));
		
		/************************** Finalize ***************************/
		printf("writing file......\n");
		int* pagerank_rank = getrank(pagerank, num_nodes);
		FILE *fp_write = fopen("./pagerank_parallel.result", "w");
		fprintf(fp_write, "time: %lu ms\n", 1000000 * (endtime.tv_sec - starttime.tv_sec) + (endtime.tv_usec - starttime.tv_usec));
		fprintf(fp_write, "node_id\tpagerank\n");
		for (i = 0; i < num_nodes; i++) {
			fprintf(fp_write, "%d\t%d\n", i, pagerank_rank[i]);
		}
		fclose(fp_write);
		printf("Finish writing file\n");
	}
	/** Slaves **/
	else {
		int data_recved;
		MPI_Status status;
		
		MPI_Recv(recv_buf, pre_count[myrank], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_DOUBLE, &data_recved);
		while (data_recved != 1) {
			for (i = 0; i < pre_count[myrank]; i++) {
				temp_pagerank[pre_mapping[myrank][i]] = recv_buf[i];
			}
			
			current = row_ptrs_partition;
			i = 0;
			k = 0;
			while (current->next != NULL) {
				if (post_mapping[myrank][k] == i) {
					double temp = 0;
					for (j = current->val; j < current->next->val; j++) {
						temp += temp_pagerank[col_inds_partition[j]] * col_prob[col_inds_partition[j]];
					}
					pagerank[k] = temp; // compresssed length
					k++;
				}
				current = current->next;
				i++;
			}
			
			
			MPI_Send(pagerank, post_count[myrank], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
			
			memset((char*)recv_buf, 0, sizeof(double) * num_nodes);
			memset((char*)temp_pagerank, 0, sizeof(double) * num_nodes);
			// memset((char*)pagerank, 0, sizeof(double) * num_nodes);
			
			MPI_Recv(recv_buf, pre_count[myrank], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_DOUBLE, &data_recved);
		}
	}
	
	MPI_Finalize();
}


