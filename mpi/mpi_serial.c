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
	int i, j, k;
	
	/****************** Read split file **********************/
	int num_nodes = 0;
	char text[32];
	char* split_info = argv[2];
	FILE *fp_split = fopen(split_info, "r");
	while (fscanf(fp_split, "%s", text) != EOF) {
		num_nodes++;
	}
	printf("num of nodes: %d\n", num_nodes);
	memset(text, 0, 32);
	
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
	
	// printf("Number of column in CRS: %d\n", col_len);
	
	/**************************** Initialize **************************/
	double* pagerank = (double*)malloc(sizeof(double) * num_nodes);
	double* temp_pagerank = (double*)malloc(sizeof(double) * num_nodes);
	
	for (i = 0; i < num_nodes; i++) {
		pagerank[i] = 1.0 / num_nodes;
	}

	struct graph_node *current = row_ptrs_head;
	double *col_prob = (double*)malloc(sizeof(double) * num_nodes);
	i = 0;
	while (current->next != NULL) {
		col_prob[i] = 1 / (double)(current->next->val - current->val);
		current = current->next;
		i++;
	}
	// printf("preparation done!\n");
	
	double last_norm = norm2(pagerank, num_nodes);
	printf("original norm: %f\n", last_norm);
	
	/************************* Calculate *************************/
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
		if (fabs(norm2(temp_pagerank, num_nodes) - last_norm) < 0.00001) {
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
	
	/************************** Finalize ***************************/
	printf("writing file......\n");
	int* pagerank_rank = getrank(pagerank, num_nodes);
	FILE *fp_write = fopen("./pagerank_serial.result", "w");
	fprintf(fp_write, "tine: %lu ms\n", 1000000 * (endtime.tv_sec - starttime.tv_sec) + (endtime.tv_usec - starttime.tv_usec));
	fprintf(fp_write, "node_id\tpagerank\n");
	for (i = 0; i < num_nodes; i++) {
		fprintf(fp_write, "%d\t%d\n", i, pagerank_rank[i]);
	}
	fclose(fp_write);
	printf("Finish writing file\n");
}


