//Copyright <2016> <Chu-Min Li & Hua Jiang>
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <limits.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <vector>
class LMCMaxClique {
#define WORD_LENGTH 100
#define TRUE 1
#define FALSE 0
#define NONE -1
#define DELIMITER 0
#define PASSIVE 0
#define ACTIVE 1
#define P_TRUE 2
#define P_FALSE 0
#define NO_REASON -3
#define CONFLICT -1978
#define MAX_NODE 100000
#define max_expand_depth 100000
#define PUSH_OPT(item, stack) stack[stack ## _fill_pointer++] = item
#define ptr(stack) stack ## _fill_pointer
#define is_neibor(i,j) matrice[i][j]

#define CUR_CLQ_SIZE Clique_Stack_fill_pointer
#define CURSOR Cursor_Stack[Cursor_Stack_fill_pointer-1]
#define MIN(a,b) a<=b?a:b
#define BIT_MAP_SIZE 4097

#define SET_EDGE(row,col) ((*(Adj_Matrix + (row)* MATRIX_ROW_WIDTH + ((col) >> 3))) |= (1 << ((col) & 7)))
#define GET_EDGE(row,col) ((*(Adj_Matrix + (row)* MATRIX_ROW_WIDTH + ((col) >> 3))) & (1 << ((col) & 7)))

#define iMatrix(i) (Adj_Matrix+(i)*MATRIX_ROW_WIDTH)
#define Matrix(i,j) ((*((i) + ((j) >> 3))) & (1 << ((j) & 7)))

    int FORMAT = 1, nodeNum, NB_NODE_O=0,  edgeNum, NB_EDGE_O=0,
        MAX_CLQ_SIZE=0, MAX_ISET_SIZE=0, INIT_CLQ_SIZE=0,MATRIX_ROW_WIDTH=0, 
        MAX_VERTEX_NO=0, K_CORE_G = 0;

#define CORE_NO Vertex_UB
    int Max_Degree = 0;
    int Node_Degree[MAX_NODE];
    char Node_State[MAX_NODE];
    int** Node_Neibors;

    int Candidate_Stack_fill_pointer = 0;
    int Candidate_Stack[MAX_NODE * 2];
    int Vertex_UB[MAX_NODE * 2];
    int Clique_Stack_fill_pointer;
    int* Clique_Stack, * MaxCLQ_Stack;
    int Cursor_Stack[max_expand_depth];
    int Cursor_Stack_fill_pointer = 0;

    int* Node_Reason;
    unsigned char* Adj_Matrix;

    int iSET_COUNT = 0;
    int* iSET_Size;
    char* iSET_State;
    char* iSET_Used;
    char* iSET_Tested;
    int* iSET_Index;
    char* iSET_Involved;
    char* Is_Tested;
    int** iSET;

    int* REASON_STACK;
    int REASON_STACK_fill_pointer = 0;
    int* CONFLICT_ISET_STACK;
    int CONFLICT_ISET_STACK_fill_pointer;
    int* ADDED_NODE_iSET;
    int* REDUCED_iSET_STACK = Node_Degree;
    int REDUCED_iSET_STACK_fill_pointer = 0;
    int* PASSIVE_iSET_STACK;
    int PASSIVE_iSET_STACK_fill_pointer = 0;
    int* FIXED_NODE_STACK;
    int FIXED_NODE_STACK_fill_pointer = 0;
    int* UNIT_STACK;
    int UNIT_STACK_fill_pointer = 0;
    int* NEW_UNIT_STACK;
    int NEW_UNIT_STACK_fill_pointer = 0;

    int Rollback_Point;
    int Branching_Point;
    //int Matrix_Reallocation = 0;
    int* Old_Name;
    int* Second_Name;
    int NB_CANDIDATE = 0, FIRST_INDEX;
    int START_MAXSAT_THD = 15;

    int Extra_Node_Stack[100000];

    int Last_Idx = 0;
    int cut_ver = 0, total_cut_ver = 0;
    int cut_inc = 0, total_cut_inc = 0;
    int cut_iset = 0, total_cut_iset = 0;
    int cut_satz = 0, total_cut_satz = 0;
    long long Branches_Nodes[6];
    int LAST_IN;
    float Dynamic_Radio = 0.70;
    int REBUILD_MATRIX = FALSE;
    int CUR_MAX_NODE;
    int Branches[1200];
    int* Init_Adj_List;
    int BLOCK_COUNT = 0;
    int* BLOCK_LIST[100];
    double READ_TIME, INIT_TIME, SEARCH_TIME;
    long long N0_0 = 0, N0_1 = 0, N1_0 = 0, N1_1 = 0, L1 = 0;
    double D0 = 0, D1 = 0, D2 = 0, D3 = 0, Dt = 0, R1 = 0;
    long long SubNB = 0;
    template<typename T>using Vec = std::vector<T>;
public:
    struct NewEdge {
        int src;
        int dst;
    };
protected:
    static int int_cmp(const void* a, const void* b) {
        return *((int*)b) - *((int*)a);
    }
    int is_adjacent(int node1, int node2) {
        int neibor, * neibors;
        neibors = Node_Neibors[node1];
        for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
            if (neibor == node2) {
                return TRUE;
            }
        }
        return FALSE;
    }

    void check_clique() {
        int i, j;
        for (i = 0; i < MAX_CLQ_SIZE; i++) {
            for (j = i + 1; j < MAX_CLQ_SIZE; j++) {
                if (MaxCLQ_Stack[i] > MaxCLQ_Stack[j]) {
                    assert(is_adjacent(MaxCLQ_Stack[i], MaxCLQ_Stack[j]));
                }
                else {
                    assert(is_adjacent(MaxCLQ_Stack[j], MaxCLQ_Stack[i]));
                }
            }
        }
    }

    void check_clique_in_result_file(int* solvers) {
        int i = 0, j, node, clique_size = 0;
        Vec<int> clique;

        clique_size = sizeof(solvers) / sizeof(solvers[0]);
        if (clique_size > 0) {
            for (i = 0; i < clique_size; i++) {
                for (j = i + 1; j < clique_size; j++) {
                    if (is_adjacent(clique[i], clique[j]) == FALSE) {
                        //printf("find non-adjacent vertices: %d %d\n", clique[i],clique[j]);
                        //printf("#FALSE:%s\n", solvers);
                        return;
                    }

                }
            }
            //printf("#TRUE:%s\n", solvers);
        }
        else {
            //printf("#NONE:%s\n", solvers);
        }
    }


    void allcoate_memory_for_adjacency_list(int nb_node, int nb_edge,
        int offset) {
        int i, block_size = 40960000, free_size = 0;
        Init_Adj_List = (int*)malloc((2 * nb_edge + nb_node) * sizeof(int));
        if (Init_Adj_List == NULL) {
            for (i = 1; i <= nodeNum; i++) {
                if (Node_Degree[i - offset] + 1 > free_size) {
                    Node_Neibors[i] = (int*)malloc(block_size * sizeof(int));
                    BLOCK_LIST[BLOCK_COUNT++] = Node_Neibors[i];
                    free_size = block_size - (Node_Degree[i - offset] + 1);
                }
                else {
                    Node_Neibors[i] = Node_Neibors[i - 1]
                        + Node_Degree[i - 1 - offset] + 1;
                    free_size = free_size - (Node_Degree[i - offset] + 1);
                }
            }
        }
        else {
            BLOCK_COUNT = 1;
            BLOCK_LIST[BLOCK_COUNT - 1] = Init_Adj_List;
            Node_Neibors[1] = Init_Adj_List;
            for (i = 2; i <= nodeNum; i++) {
                Node_Neibors[i] = Node_Neibors[i - 1] + Node_Degree[i - 1 - offset]
                    + 1;
            }
        }
    }


    int sort_by_degeneracy_ordering() {
        int* degree_counter, * where;
        int p1, i, node, neibor, * neibors, t, j, h, k;
        int cur_degree = 0;
        INIT_CLQ_SIZE = 0;
        ////printf("I computing initial degeneracy ordering...\n");
        //printf("the initl MAX_VERTEX_NO:%d******************\n", MAX_VERTEX_NO);
        where = Candidate_Stack + nodeNum + 1;
        degree_counter = Vertex_UB + nodeNum + 1;
        memset(degree_counter, 0, (Max_Degree + 1) * sizeof(int));

        for (node = 1; node <= nodeNum; node++) {
            degree_counter[Node_Degree[node]]++;
        }
        j = 0;
        for (i = 0; i <= Max_Degree; i++) {
            k = degree_counter[i];
            degree_counter[i] = j;
            j += k;
        }

        for (node = 1; node <= nodeNum; node++) {
            Candidate_Stack[t = degree_counter[Node_Degree[node]]++] = node;
            where[node] = t;
        }
        for (i = Max_Degree; i > 0; i--) {
            degree_counter[i] = degree_counter[i - 1];
        }
        degree_counter[0] = 0;

        Candidate_Stack[nodeNum] = DELIMITER;
        ptr(Candidate_Stack) = nodeNum + 1;

        p1 = 0;
        cur_degree = Node_Degree[Candidate_Stack[p1]];
        while (p1 < nodeNum) {
            node = Candidate_Stack[p1];

            if (Node_Degree[node] > cur_degree)
                cur_degree = Node_Degree[node];
            CORE_NO[p1] = cur_degree;

            if (cur_degree > K_CORE_G)
                K_CORE_G = cur_degree;

            if (p1 < nodeNum - 1
                && Node_Degree[node] == Node_Degree[Candidate_Stack[p1 + 1]]) {
                degree_counter[Node_Degree[node]] = p1 + 1;
            }
            if (Node_Degree[node] > MAX_VERTEX_NO)
                MAX_VERTEX_NO = Node_Degree[node];

            if (Node_Degree[node] == nodeNum - p1 - 1) {
                INIT_CLQ_SIZE = nodeNum - p1;
                //HEUR_CLQ_SIZE = INIT_CLQ_SIZE;
                //printf("I initial clique is %d\n", INIT_CLQ_SIZE);
                //printf("I maxcore number is %d\n", K_CORE_G);
                MaxCLQ_Stack = (int*)malloc((K_CORE_G + 2) * sizeof(int));
                Clique_Stack = (int*)malloc((K_CORE_G + 2) * sizeof(int));
                memcpy(MaxCLQ_Stack, Candidate_Stack + p1,
                    INIT_CLQ_SIZE * sizeof(int));
                for (i = p1 + 1; i < nodeNum; i++)
                    CORE_NO[i] = cur_degree;
                break;
            }

            neibors = Node_Neibors[node];
            for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                if (where[neibor] > p1) {
                    t = where[neibor];
                    h = degree_counter[Node_Degree[neibor]];

                    k = Candidate_Stack[h];

                    Candidate_Stack[h] = neibor;
                    where[neibor] = h;

                    Candidate_Stack[t] = k;
                    where[k] = t;

                    degree_counter[Node_Degree[neibor]]++;

                    Node_Degree[neibor]--;
                    if (Node_Degree[neibor]
                        != Node_Degree[Candidate_Stack[h - 1]]) {
                        degree_counter[Node_Degree[neibor]] = h;
                    }
                }
            }
            p1++;
        }
        //std::cout << "=============init:" << MAX_VERTEX_NO + 1 << "," << INIT_CLQ_SIZE << std::endl;
        if (MAX_VERTEX_NO + 1 == INIT_CLQ_SIZE) {
            MAX_CLQ_SIZE = INIT_CLQ_SIZE;
            //printf("I find the maximum clique in initial phase!\n");
            return TRUE;
        }
        else {
            return FALSE;
        }
    }

    int re_number_adj(int node) {
        int i, k, * neibors, * saved_neibors = NULL, neibor, one_neibor;
        for (i = 0; i < iSET_COUNT - 1; i++) {
            neibors = iSET[i];
            one_neibor = NONE;
            for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                if (is_adjacent(node, neibor) > 0) {
                    if (one_neibor == NONE) {
                        one_neibor = neibor;
                        saved_neibors = neibors;
                    }
                    else if (one_neibor != NONE) {
                        break;
                    }
                }
            }
            if (one_neibor == NONE) {
                iSET[i][iSET_Size[i]] = node;
                iSET_Size[i]++;
                iSET[i][iSET_Size[i]] = NONE;
                return TRUE;
            }
            if (neibor == NONE) {
                for (k = i + 1; k < iSET_COUNT; k++) {
                    neibors = iSET[k];
                    for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                        if (one_neibor > neibor) {
                            if (is_adjacent(one_neibor, neibor) == TRUE)
                                break;
                        }
                        else {
                            if (is_adjacent(neibor, one_neibor) == TRUE)
                                break;
                        }
                    }
                    if (neibor == NONE) {
                        iSET[k][iSET_Size[k]] = one_neibor;
                        iSET_Size[k]++;
                        iSET[k][iSET_Size[k]] = NONE;
                        *saved_neibors = node;
                        return TRUE;
                    }
                }
            }
        }
        return FALSE;
    }

    int re_number(int node) {
        int i, k, * neibors, * saved_neibors = NULL, neibor, one_neibor;
        unsigned char* adj_list1 = iMatrix(node), * adj_list2;
        for (i = 0; i < iSET_COUNT - 1; i++) {
            neibors = iSET[i];
            one_neibor = NONE;
            for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                if (Matrix(adj_list1, neibor) > 0) {
                    if (one_neibor == NONE) {
                        one_neibor = neibor;
                        saved_neibors = neibors;
                    }
                    else if (one_neibor != NONE) {
                        break;
                    }
                }
            }
            if (one_neibor == NONE) {
                iSET[i][iSET_Size[i]] = node;
                iSET_Size[i]++;
                iSET[i][iSET_Size[i]] = NONE;
                return TRUE;
            }
            if (neibor == NONE) {
                adj_list2 = iMatrix(one_neibor);
                for (k = i + 1; k < iSET_COUNT; k++) {
                    neibors = iSET[k];
                    for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                        if (Matrix(adj_list2, neibor) > 0)
                            break;
                    }
                    if (neibor == NONE) {
                        iSET[k][iSET_Size[k]] = one_neibor;
                        iSET_Size[k]++;
                        iSET[k][iSET_Size[k]] = NONE;
                        iSET_Index[one_neibor] = k;
                        *saved_neibors = node;
                        iSET_Index[node] = i;
                        return TRUE;
                    }
                }
            }
        }
        return FALSE;
    }

    int addIntoIsetTomitaBis_adj(int node) {
        int j, * current_set, iset_node;
        for (j = 0; j < iSET_COUNT; j++) {
            current_set = iSET[j];
            for (iset_node = *current_set; iset_node != NONE; iset_node =
                *(++current_set)) {
                //assert(node > iset_node);
                if (is_adjacent(node, iset_node) == TRUE)
                    break;
            }
            if (iset_node == NONE) {
                iSET_Size[j]++;
                *(current_set) = node;
                *(++current_set) = NONE;
                return TRUE;
            }
        }
        if (iSET_COUNT < MAX_ISET_SIZE) {
            iSET_Size[j] = 1;
            iSET[j][0] = node;
            iSET[j][1] = NONE;
            iSET_COUNT++;
            return TRUE;
        }
        else {
            return FALSE;
        }
    }

    int addIntoIsetTomitaBis(int node) {
        int j, * current_set, iset_node;
        unsigned char* adj_list = iMatrix(node);
        for (j = 0; j < iSET_COUNT; j++) {
            current_set = iSET[j];
            for (iset_node = *current_set; iset_node != NONE; iset_node =
                *(++current_set)) {
                if (Matrix(adj_list, iset_node) > 0)
                    break;
            }
            if (iset_node == NONE) {
                iSET_Size[j]++;
                *(current_set) = node;
                *(++current_set) = NONE;
                iSET_Index[node] = j;
                return TRUE;
            }
        }
        if (iSET_COUNT < MAX_ISET_SIZE) {
            iSET_Size[j] = 1;
            iSET[j][0] = node;
            iSET[j][1] = NONE;
            iSET_Index[node] = j;
            iSET_COUNT++;
            return TRUE;
        }
        else {
            return FALSE;
        }
    }

    int cut_by_iset_less_vertices() {
        int i = ptr(Candidate_Stack) - 2, node;
        iSET_COUNT = 0;
        MAX_ISET_SIZE = MAX_CLQ_SIZE - CUR_CLQ_SIZE - 1;
        for (node = Candidate_Stack[i]; node != DELIMITER; node =
            Candidate_Stack[--i]) {
            if (addIntoIsetTomitaBis_adj(node) == FALSE
                && re_number_adj(node) == FALSE) {
                return FALSE;
            }
        }
        return TRUE;
    }

    int cut_by_iset_last_renumber() {
        int i = ptr(Candidate_Stack) - 2, node;
        LAST_IN = INT_MAX;
        FIRST_INDEX = NONE;
        iSET_COUNT = 0;
        MAX_ISET_SIZE = MAX_CLQ_SIZE - CUR_CLQ_SIZE - 1;
        for (node = Candidate_Stack[i]; node != DELIMITER; node =
            Candidate_Stack[--i]) {
            if (addIntoIsetTomitaBis(node) == FALSE && re_number(node) == FALSE) {
                if (FIRST_INDEX == NONE)
                    FIRST_INDEX = i;
                Candidate_Stack[i] = -node;
            }
            else {
                LAST_IN = i;
            }
        }
        if (FIRST_INDEX == NONE) {
            cut_iset++;
            Vertex_UB[CURSOR] = iSET_COUNT + 1;
            return TRUE;
        }
        else {
            Branching_Point = FIRST_INDEX + 1;
            if (MAX_CLQ_SIZE < START_MAXSAT_THD) {
                i = FIRST_INDEX;
                for (node = Candidate_Stack[i]; node != DELIMITER; node =
                    Candidate_Stack[--i]) {
                    if (node < 0 && i > LAST_IN) {
                        Vertex_UB[i] = MAX_ISET_SIZE + 1;
                    }
                }
            }
            return FALSE;
        }
    }


#define assign_node(node, value, reason) {\
	Node_State[node] = value;\
	Node_Reason[node] = reason;\
	PUSH_OPT(node, FIXED_NODE_STACK);\
}
    int fix_newNode_for_iset(int fix_node, int fix_iset) {
        int idx, iset_idx;
        iSET_State[fix_iset] = PASSIVE;
        PUSH_OPT(fix_iset, PASSIVE_iSET_STACK);
        assign_node(fix_node, P_FALSE, fix_iset);
        //iSET_Potential[fix_iset] -= Node_Potential[fix_node];
        idx = ADDED_NODE_iSET[fix_node];
        while ((iset_idx = CONFLICT_ISET_STACK[idx++]) != NONE) {
            if (iSET_State[iset_idx] == ACTIVE) {
                iSET_Size[iset_idx]--;
                //iSET_ADD_Size[iset_idx]--;
                PUSH_OPT(iset_idx, REDUCED_iSET_STACK);
                //iSET_Potential[iset_idx] -= Node_Potential[fix_node];
                if (iSET_Size[iset_idx] == 1)
                    PUSH_OPT(iset_idx, NEW_UNIT_STACK);
                else if (iSET_Size[iset_idx] == 0)
                    return iset_idx;
            }
        }
        return NONE;
    }

    int fix_oldNode_for_iset(int fix_node, int fix_iset) {
        int i, node, iset_idx;
        unsigned char* adj_list;
        assign_node(fix_node, P_TRUE, fix_iset);
        iSET_State[fix_iset] = PASSIVE;
        PUSH_OPT(fix_iset, PASSIVE_iSET_STACK);

        adj_list = iMatrix(fix_node);
        i = ptr(Candidate_Stack) - 2;
        for (node = Candidate_Stack[i]; node != DELIMITER; node =
            Candidate_Stack[--i]) {
            if (node
            > 0 && Matrix(adj_list, node) == FALSE && Node_State[node] == ACTIVE) {
                assign_node(node, P_FALSE, fix_iset);
                iset_idx = iSET_Index[node];
                if (iSET_State[iset_idx] == ACTIVE) {
                    iSET_Size[iset_idx]--;
                    PUSH_OPT(iset_idx, REDUCED_iSET_STACK);
                    if (iSET_Size[iset_idx] == 1) {
                        PUSH_OPT(iset_idx, NEW_UNIT_STACK);
                    }
                    else if (iSET_Size[iset_idx] == 0) {
                        return iset_idx;
                    }
                }
            }
        }

        return NONE;
    }

#define fix_node(node,iset) ((node > nodeNum)? fix_newNode_for_iset(node, iset):fix_oldNode_for_iset(node, iset))

    int fix_node_iset(int fix_iset) {
        int fix_node, * nodes = iSET[fix_iset];
        for (fix_node = *(nodes); fix_node != NONE; fix_node = *(++nodes)) {
            if (Node_State[fix_node] == ACTIVE) {
                if (fix_node > MAX_VERTEX_NO)
                    return fix_newNode_for_iset(fix_node, fix_iset);
                else
                    return fix_oldNode_for_iset(fix_node, fix_iset);
            }
        }
        nodes = iSET[fix_iset];
        for (fix_node = *(nodes); fix_node != NONE; fix_node = *(++nodes)) {
            //printf("iset=%d,node=%d,active=%d\n", fix_iset, fix_node,Node_State[fix_node]);
        }
        //printf("error in fix_node_iset\n");
        //printf("iSET COUNT=%d\n", iSET_COUNT);
        //print_all_iset();
        exit(0);
    }

    int unit_iset_process() {
        int i = 0, j = 0, iset_idx, empty_iset;
        for (i = 0; i < ptr(UNIT_STACK); i++) {
            iset_idx = UNIT_STACK[i];
            if (iSET_State[iset_idx] == ACTIVE && iSET_Size[iset_idx] == 1) {
                ptr(NEW_UNIT_STACK) = 0;
                if ((empty_iset = fix_node_iset(iset_idx)) > NONE) {
                    return empty_iset;
                }
                else {
                    for (j = 0; j < ptr(NEW_UNIT_STACK); j++) {
                        iset_idx = NEW_UNIT_STACK[j];
                        if (iSET_State[iset_idx] == ACTIVE) {
                            if ((empty_iset = fix_node_iset(iset_idx)) > NONE)
                                return empty_iset;
                        }
                    }
                }
            }
        }
        ptr(NEW_UNIT_STACK) = 0;
        return NONE;
    }

    int unit_iset_process_used_first() {
        int j, iset, iset_start = 0, used_iset_start = 0, my_iset;
        do {
            for (j = used_iset_start; j < ptr(NEW_UNIT_STACK); j++) {
                iset = NEW_UNIT_STACK[j];
                if (iSET_State[iset] == ACTIVE && iSET_Used[iset] == TRUE)
                    if ((my_iset = fix_node_iset(iset)) != NONE)
                        return my_iset;
            }
            used_iset_start = j;
            for (j = iset_start; j < ptr(NEW_UNIT_STACK); j++) {
                iset = NEW_UNIT_STACK[j];
                if (iSET_State[iset] == ACTIVE) {
                    if ((my_iset = fix_node_iset(iset)) != NONE)
                        return my_iset;
                    iset_start = j + 1;
                    break;
                }
            }
        } while (j < ptr(NEW_UNIT_STACK));
        return NONE;
    }


    void identify_conflict_sets(int iset_idx) {
        int i, reason_start = ptr(REASON_STACK), iset, * nodes, node, reason_iset;
        PUSH_OPT(iset_idx, REASON_STACK);
        iSET_Involved[iset_idx] = TRUE;
        for (i = reason_start; i < ptr(REASON_STACK); i++) {
            iset = REASON_STACK[i];
            nodes = iSET[iset];
            for (node = *nodes; node != NONE; node = *(++nodes))
                if (Node_State[node] == P_FALSE && Node_Reason[node] != NO_REASON
                    && iSET_Involved[Node_Reason[node]] == FALSE) {
                    reason_iset = Node_Reason[node];
                    PUSH_OPT(reason_iset, REASON_STACK);
                    Node_Reason[node] = NO_REASON;
                    iSET_Involved[reason_iset] = TRUE;
                }
        }
        for (i = reason_start; i < ptr(REASON_STACK); i++) {
            iSET_Involved[REASON_STACK[i]] = FALSE;
            iSET_Used[REASON_STACK[i]] = TRUE;
        }
    }
    void enlarge_conflict_sets(int & ADDED_NODE) {
        int i, iset;
        Node_State[ADDED_NODE] = ACTIVE;
        Node_Reason[ADDED_NODE] = NO_REASON;
        //node_match_state[ADDED_NODE] = FALSE;
        ADDED_NODE_iSET[ADDED_NODE] = ptr(CONFLICT_ISET_STACK);
        for (i = 0; i < ptr(REASON_STACK); i++) {
            iset = REASON_STACK[i];
            if (iSET_Involved[iset] == FALSE) {
                iSET_Involved[iset] = TRUE;
                iSET[iset][iSET_Size[iset]++] = ADDED_NODE;
                iSET[iset][iSET_Size[iset]] = NONE;
                //assert(ptr(CONFLICT_ISET_STACK)<3*tab_node_size);
                PUSH_OPT(iset, CONFLICT_ISET_STACK);
            }
        }
        PUSH_OPT(NONE, CONFLICT_ISET_STACK);
        i = ADDED_NODE_iSET[ADDED_NODE];
        //	while ((iset = CONFLICT_ISET_STACK[i++]) != NONE) {
        //		iSET_Potential[iset] += Node_Potential[ADDED_NODE];
        //	}
        for (i = 0; i < ptr(REASON_STACK); i++) {
            iSET_Involved[REASON_STACK[i]] = FALSE;
            iSET_Used[REASON_STACK[i]] = FALSE;
        }
        ptr(REASON_STACK) = 0;
        ADDED_NODE++;
    }

    void rollback_context_for_maxsatz(int start_fixed, int start_passive,
        int start_reduced) {
        int i, node;
        for (i = start_fixed; i < ptr(FIXED_NODE_STACK); i++) {
            node = FIXED_NODE_STACK[i];
            Node_State[node] = ACTIVE;
            Node_Reason[node] = NO_REASON;
        }
        ptr(FIXED_NODE_STACK) = start_fixed;
        for (i = start_passive; i < ptr(PASSIVE_iSET_STACK); i++) {
            iSET_State[PASSIVE_iSET_STACK[i]] = ACTIVE;
        }
        ptr(PASSIVE_iSET_STACK) = start_passive;
        for (i = start_reduced; i < ptr(REDUCED_iSET_STACK); i++) {
            iSET_Size[REDUCED_iSET_STACK[i]]++;
        }
        ptr(REDUCED_iSET_STACK) = start_reduced;
        ptr(NEW_UNIT_STACK) = 0;
    }

    void reset_context_for_maxsatz() {
        int i, node;
        for (i = 0; i < ptr(FIXED_NODE_STACK); i++) {
            node = FIXED_NODE_STACK[i];
            Node_State[node] = ACTIVE;
            Node_Reason[node] = NO_REASON;
        }
        ptr(FIXED_NODE_STACK) = 0;
        for (i = 0; i < ptr(PASSIVE_iSET_STACK); i++) {
            iSET_State[PASSIVE_iSET_STACK[i]] = ACTIVE;
        }
        ptr(PASSIVE_iSET_STACK) = 0;
        for (i = 0; i < ptr(REDUCED_iSET_STACK); i++) {
            iSET_Size[REDUCED_iSET_STACK[i]]++;
        }
        ptr(REDUCED_iSET_STACK) = 0;
        ptr(NEW_UNIT_STACK) = 0;
    }

    int further_test_reduced_iset(int start) {
        int i, chosen_iset, empty_iset, node, * nodes;
        int saved_fixed = ptr(FIXED_NODE_STACK);
        int start_passive = ptr(PASSIVE_iSET_STACK);
        int start_reduced = ptr(REDUCED_iSET_STACK);
        for (i = start; i < ptr(REDUCED_iSET_STACK); i++)
            iSET_Tested[REDUCED_iSET_STACK[i]] = FALSE;
        for (i = start; i < ptr(REDUCED_iSET_STACK); i++) {
            chosen_iset = REDUCED_iSET_STACK[i];
            if (iSET_Tested[chosen_iset] == FALSE
                && iSET_State[chosen_iset] == ACTIVE
                && iSET_Size[chosen_iset] == 2) {
                iSET_Tested[chosen_iset] = TRUE;
                nodes = iSET[chosen_iset];
                for (node = *nodes; node != NONE; node = *(++nodes)) {
                    if (Node_State[node] == ACTIVE && node > MAX_VERTEX_NO)
                        break;
                }
                if (node == NONE) {
                    nodes = iSET[chosen_iset];
                    for (node = *nodes; node != NONE; node = *(++nodes)) {
                        if (Node_State[node] == ACTIVE) {
                            ptr(NEW_UNIT_STACK) = 0;
                            if ((empty_iset = fix_oldNode_for_iset(node,
                                chosen_iset)) != NONE || (empty_iset =
                                    unit_iset_process_used_first()) != NONE) {
                                iSET_Involved[chosen_iset] = TRUE;
                                identify_conflict_sets(empty_iset);
                                iSET_Involved[chosen_iset] = FALSE;
                                rollback_context_for_maxsatz(saved_fixed,
                                    start_passive, start_reduced);
                            }
                            else {
                                rollback_context_for_maxsatz(saved_fixed,
                                    start_passive, start_reduced);
                                break;
                            }
                        }
                    }
                    if (node == NONE)
                        return chosen_iset;
                }
            }
        }
        return NONE;
    }

    int fix_anyNode_for_iset(int fix_node, int fix_iset) {
        if (fix_node > MAX_VERTEX_NO)
            return fix_newNode_for_iset(fix_node, fix_iset);
        else
            return fix_oldNode_for_iset(fix_node, fix_iset);
    }

    int inc_maxsatz_lookahead_by_fl2() {
        int i, j, empty_iset, iset_idx, * nodes, node;
        int sn = FIXED_NODE_STACK_fill_pointer;
        int sp = PASSIVE_iSET_STACK_fill_pointer;
        int sr = REDUCED_iSET_STACK_fill_pointer;
        int rs = REASON_STACK_fill_pointer;
        for (i = 0; i < sr; i++) {
            Is_Tested[REDUCED_iSET_STACK[i]] = FALSE;
        }
        for (i = 0; i < sr; i++) {
            iset_idx = REDUCED_iSET_STACK[i];
            if (Is_Tested[iset_idx] == FALSE && iSET_State[iset_idx] == ACTIVE
                && iSET_Size[iset_idx] == 2) {
                nodes = iSET[iset_idx];
                ptr(REASON_STACK) = rs;
                Is_Tested[iset_idx] = TRUE;
                for (node = *nodes; node != NONE; node = *(++nodes)) {
                    if (Node_State[node] == ACTIVE) {
                        ptr(NEW_UNIT_STACK) = 0;
                        if ((empty_iset = fix_anyNode_for_iset(node, iset_idx))
                            != NONE || (empty_iset =
                                unit_iset_process_used_first()) != NONE
                            || (*(nodes + 1) == NONE && (empty_iset =
                                further_test_reduced_iset(sr)) != NONE)) {
                            identify_conflict_sets(empty_iset);
                            rollback_context_for_maxsatz(sn, sp, sr);
                        }
                        else {
                            rollback_context_for_maxsatz(sn, sp, sr);
                            break;
                        }
                    }
                }
                if (node == NONE) {
                    //reset_context_for_maxsatz();
                    return TRUE;
                }
                else {
                    for (j = rs; j < ptr(REASON_STACK); j++) {
                        iSET_Involved[REASON_STACK[j]] = FALSE;
                        iSET_Used[REASON_STACK[j]] = FALSE;
                    }
                    ptr(REASON_STACK) = rs;
                }
            }
        }
        //reset_context_for_maxsatz();
        return FALSE;
    }

    int inc_maxsatz_on_last_iset(int &ADDED_NODE) {
        int empty_iset, node, * nodes, iset_idx = iSET_COUNT - 1;
        ptr(REASON_STACK) = 0;
        ptr(NEW_UNIT_STACK) = 0;
        ptr(FIXED_NODE_STACK) = 0;
        ptr(PASSIVE_iSET_STACK) = 0;
        ptr(REDUCED_iSET_STACK) = 0;
        nodes = iSET[iset_idx];
        for (node = *nodes; node != NONE; node = *(++nodes)) {
            if (Node_State[node] == ACTIVE) {
                ptr(NEW_UNIT_STACK) = 0;
                if ((empty_iset = fix_oldNode_for_iset(node, iset_idx)) != NONE
                    || (empty_iset = unit_iset_process_used_first()) != NONE
                    || (empty_iset = unit_iset_process()) != NONE) {
                    identify_conflict_sets(empty_iset);
                    reset_context_for_maxsatz();
                }
                else if (inc_maxsatz_lookahead_by_fl2() == TRUE) {
                    reset_context_for_maxsatz();
                }
                else {
                    reset_context_for_maxsatz();
                    break;
                }
            }
        }
        if (node == NONE) {
            enlarge_conflict_sets(ADDED_NODE);
        }
        return node;
    }

    int open_new_iset_old(int i) {
        int* current_set, node, iset_node, idx = 0;
        unsigned char* adj_list;
        iSET_Size[iSET_COUNT] = 0;
        iSET[iSET_COUNT][0] = NONE;
        iSET_Used[iSET_COUNT] = FALSE;
        iSET_State[iSET_COUNT] = ACTIVE;

        while ((node = Candidate_Stack[i]) != DELIMITER) {
            if (node < 0) {
                node = -node;
                adj_list = iMatrix(node);
                current_set = iSET[iSET_COUNT];
                for (iset_node = *current_set; iset_node != NONE; iset_node =
                    *(++current_set)) {
                    if (Matrix(adj_list, iset_node) > 0)
                        break;
                }
                if (iset_node == NONE) {
                    iSET_Size[iSET_COUNT]++;
                    *(current_set) = node;
                    *(++current_set) = NONE;
                    iSET_Index[node] = iSET_COUNT;
                    Node_State[node] = ACTIVE;
                    Node_Reason[node] = NO_REASON;
                    Candidate_Stack[i] = node;
                    Extra_Node_Stack[idx++] = node;
                    Extra_Node_Stack[idx++] = i;
                }
                else {
                    break;
                }
            }
            i--;
        }
        if (iSET_Size[iSET_COUNT] == 1) {
            PUSH_OPT(UNIT_STACK[0], UNIT_STACK);
            UNIT_STACK[0] = iSET_COUNT;
        }
        iSET_COUNT++;
        //assert(iSET_COUNT <= MAX_VERTEX_NO);
        Extra_Node_Stack[idx] = NONE;
        return i;
    }


    int simple_further_test_node(int start) {
        int my_iset, saved_node_stack_fill_pointer,
            saved_passive_iset_stack_fill_pointer,
            saved_reduced_iset_stack_fill_pointer, chosen_iset, node, * nodes, i,
            j;
        int my_saved_node_stack_fill_pointer,
            my_saved_passive_iset_stack_fill_pointer,
            my_saved_reduced_iset_stack_fill_pointer, conflict = FALSE;

        saved_node_stack_fill_pointer = FIXED_NODE_STACK_fill_pointer;
        saved_passive_iset_stack_fill_pointer = PASSIVE_iSET_STACK_fill_pointer;
        saved_reduced_iset_stack_fill_pointer = REDUCED_iSET_STACK_fill_pointer;
        my_saved_node_stack_fill_pointer = FIXED_NODE_STACK_fill_pointer;
        my_saved_passive_iset_stack_fill_pointer = PASSIVE_iSET_STACK_fill_pointer;
        my_saved_reduced_iset_stack_fill_pointer = REDUCED_iSET_STACK_fill_pointer;

        for (i = start; i < ptr(REDUCED_iSET_STACK); i++)
            iSET_Tested[REDUCED_iSET_STACK[i]] = FALSE;
        for (i = start; i < ptr(REDUCED_iSET_STACK); i++) {
            chosen_iset = REDUCED_iSET_STACK[i];
            if (iSET_State[chosen_iset] == ACTIVE
                && iSET_Tested[chosen_iset] == FALSE
                && iSET_Size[chosen_iset] <= 2) {
                nodes = iSET[chosen_iset];
                iSET_Tested[chosen_iset] = TRUE;
                for (node = *nodes; node != NONE; node = *(++nodes)) {
                    if (node <= MAX_VERTEX_NO && Node_State[node] == ACTIVE) {
                        ptr(NEW_UNIT_STACK) = 0;
                        my_iset = fix_oldNode_for_iset(node, chosen_iset);
                        if (my_iset == NONE)
                            my_iset = unit_iset_process_used_first();
                        rollback_context_for_maxsatz(
                            my_saved_node_stack_fill_pointer,
                            my_saved_passive_iset_stack_fill_pointer,
                            my_saved_reduced_iset_stack_fill_pointer);
                        if (my_iset != NONE) {
                            assign_node(node, P_FALSE, NO_REASON);
                            iSET_Size[chosen_iset]--;
                            PUSH_OPT(chosen_iset, REDUCED_iSET_STACK);
                            if (iSET_Size[chosen_iset] == 1) {
                                ptr(NEW_UNIT_STACK) = 0;
                                PUSH_OPT(chosen_iset, NEW_UNIT_STACK);
                                if (unit_iset_process_used_first() != NONE) {
                                    conflict = TRUE;
                                    break;
                                }
                                for (j = my_saved_reduced_iset_stack_fill_pointer;
                                    j < ptr(REDUCED_iSET_STACK); j++)
                                    iSET_Tested[REDUCED_iSET_STACK[j]] = FALSE;
                                my_saved_node_stack_fill_pointer =
                                    FIXED_NODE_STACK_fill_pointer;
                                my_saved_passive_iset_stack_fill_pointer =
                                    PASSIVE_iSET_STACK_fill_pointer;
                                my_saved_reduced_iset_stack_fill_pointer =
                                    REDUCED_iSET_STACK_fill_pointer;
                            }
                        }
                    }
                }
                if (conflict == TRUE)
                    break;
            }
        }
        rollback_context_for_maxsatz(saved_node_stack_fill_pointer,
            saved_passive_iset_stack_fill_pointer,
            saved_reduced_iset_stack_fill_pointer);
        if (conflict == TRUE)
            return chosen_iset;
        else
            return NONE;
    }

    int test_node_for_failed_nodes(int node, int iset) {
        int my_iset, saved_node_stack_fill_pointer,
            saved_passive_iset_stack_fill_pointer,
            saved_reduced_iset_stack_fill_pointer;

        saved_node_stack_fill_pointer = FIXED_NODE_STACK_fill_pointer;
        saved_passive_iset_stack_fill_pointer = PASSIVE_iSET_STACK_fill_pointer;
        saved_reduced_iset_stack_fill_pointer = REDUCED_iSET_STACK_fill_pointer;
        ptr(NEW_UNIT_STACK) = 0;
        if ((my_iset = fix_oldNode_for_iset(node, iset)) == NONE) {
            if ((my_iset = unit_iset_process_used_first()) == NONE) {
                my_iset = simple_further_test_node(
                    saved_reduced_iset_stack_fill_pointer);
            }
        }
        rollback_context_for_maxsatz(saved_node_stack_fill_pointer,
            saved_passive_iset_stack_fill_pointer,
            saved_reduced_iset_stack_fill_pointer);
        return my_iset;
    }
    int test_by_eliminate_failed_nodes() {
        int node, my_iset, * nodes, conflict, false_flag;
        do {
            false_flag = 0;
            for (my_iset = iSET_COUNT - 1; my_iset >= 0; my_iset--) {
                if (iSET_State[my_iset] == ACTIVE) {
                    nodes = iSET[my_iset];
                    conflict = FALSE;
                    ptr(NEW_UNIT_STACK) = 0;
                    for (node = *nodes; node != NONE; node = *(++nodes)) {
                        if (node <= MAX_VERTEX_NO && Node_State[node] == ACTIVE
                            && test_node_for_failed_nodes(node, my_iset)
                            != NONE) {
                            ptr(NEW_UNIT_STACK) = 0;
                            assign_node(node, P_FALSE, NO_REASON);
                            false_flag++;
                            iSET_Size[my_iset]--;
                            PUSH_OPT(my_iset, REDUCED_iSET_STACK);
                            if (iSET_Size[my_iset] == 1) {
                                PUSH_OPT(my_iset, NEW_UNIT_STACK);
                                break;
                            }
                            else if (iSET_Size[my_iset] == 0) {
                                conflict = TRUE;
                                break;
                            }
                        }
                    }
                    if (conflict == TRUE)
                        break;
                    else if (ptr(NEW_UNIT_STACK)
                                    > 0 && unit_iset_process_used_first() != NONE) {
                        conflict = TRUE;
                        break;
                    }
                }
            }
        } while (false_flag > 1 && conflict == FALSE);
        reset_context_for_maxsatz();
        return conflict;
    }


    int cut_by_inc_maxsat_eliminate_first() {
        int i, j, k = 0, node;
        int ADDED_NODE;
        ADDED_NODE = MAX_VERTEX_NO + 1;
        ptr(CONFLICT_ISET_STACK) = 0;
        ptr(UNIT_STACK) = 0;
        cut_satz++;

        i = ptr(Candidate_Stack) - 2;
        for (node = Candidate_Stack[i]; node != DELIMITER; node =
            Candidate_Stack[--i]) {
            if (node > 0)
                Node_State[node] = ACTIVE;
            else
                Node_State[-node] = PASSIVE;
        }

        for (i = 0; i < iSET_COUNT; i++) {
            if (iSET_Size[i] == 1)
                PUSH_OPT(i, UNIT_STACK);
        }

        while ((node = Candidate_Stack[FIRST_INDEX]) != DELIMITER) {
            j = open_new_iset_old(FIRST_INDEX);
            if (Candidate_Stack[j] == DELIMITER
                && test_by_eliminate_failed_nodes() == TRUE) {
                FIRST_INDEX = j;
                break;
            }
            else if ((node = inc_maxsatz_on_last_iset(ADDED_NODE)) == NONE) {
                FIRST_INDEX = j;
            }
            else {
                for (k = 0; Extra_Node_Stack[k] != node; k += 2)
                    ;
                FIRST_INDEX = Extra_Node_Stack[k + 1];
                Branching_Point = FIRST_INDEX + 1;
                //SoMC
    //			for (k = FIRST_INDEX; Candidate_Stack[k] != DELIMITER; k--) {
    //				if (Candidate_Stack[k] > 0)
    //					Candidate_Stack[k] = -Candidate_Stack[k];
    //			}
                //DoMC
                for (; Extra_Node_Stack[k] != NONE; k += 2)
                    Candidate_Stack[Extra_Node_Stack[k + 1]] = -Extra_Node_Stack[k];

                for (k = FIRST_INDEX; Candidate_Stack[k] != DELIMITER; k--) {
                    if (Candidate_Stack[k] < 0 && k > LAST_IN)
                        Vertex_UB[k] = MAX_ISET_SIZE + 1;
                }
                cut_satz--;
                break;
            }
        }

        i = ptr(Candidate_Stack) - 2;
        while ((node = Candidate_Stack[i--]) != DELIMITER) {
            if (node > 0)
                Node_State[node] = PASSIVE;
            else
                Node_State[-node] = PASSIVE;
        }

        for (node = MAX_VERTEX_NO; node <= ADDED_NODE; node++) {
            Node_State[node] = PASSIVE;
        }

        if (Candidate_Stack[FIRST_INDEX] == DELIMITER) {
            Vertex_UB[CURSOR] = MAX_CLQ_SIZE - CUR_CLQ_SIZE;
            return TRUE;
        }
        return FALSE;
    }

    float compute_subgraph_density(int start) {
        int i = start, node, neibor, * neibors;
        int nb_node = 0, nb_edge = 0;

        for (node = Candidate_Stack[i = start]; node != DELIMITER; node =
            Candidate_Stack[++i]) {
            Node_State[node] = 99;
            nb_node++;
        }
        for (node = Candidate_Stack[i = start]; node != DELIMITER; node =
            Candidate_Stack[++i]) {
            neibors = Node_Neibors[node];
            for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                if (Node_State[neibor] == 99) {
                    nb_edge++;
                }
            }
        }

        for (node = Candidate_Stack[i = start]; node != DELIMITER; node =
            Candidate_Stack[++i]) {
            Node_State[node] = PASSIVE;
        }
        return ((float)nb_edge * 2 / nb_node / (nb_node - 1));;
    }

    int compute_subgraph_degree(int start) {
        int i = start, j = 0, node, neibor, * neibors;
        Max_Degree = 0;
        int nb_node = 0, nb_edge = 0;

        for (node = Candidate_Stack[i = start]; node != DELIMITER; node =
            Candidate_Stack[++i]) {
            Node_State[node] = ACTIVE;
            Node_Degree[node] = 0;
            iSET_Index[node] = j;
            iSET_Size[j++] = 0;
            nb_node++;
        }
        for (node = Candidate_Stack[i = start]; node != DELIMITER; node =
            Candidate_Stack[++i]) {
            neibors = Node_Neibors[node];
            for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                if (Node_State[neibor] == ACTIVE) {
                    nb_edge++;
                    iSET[iSET_Index[node]][Node_Degree[node]++] = neibor;
                    iSET[iSET_Index[neibor]][Node_Degree[neibor]++] = node;
                }
            }
        }

        for (node = Candidate_Stack[i = start]; node != DELIMITER; node =
            Candidate_Stack[++i]) {
            iSET[iSET_Index[node]][Node_Degree[node]] = NONE;
            Node_State[node] = PASSIVE;
            if (Node_Degree[node] > Max_Degree)
                Max_Degree = Node_Degree[node];
        }
        Dt = ((float)nb_edge * 2 / nb_node / (nb_node - 1));
        return FALSE;
    }
    int BRANCHING_COUNT = 0;
    void allocate_memory_for_maxsat() {
        Node_Reason = (int*)malloc((MAX_VERTEX_NO + 1) * 10 * sizeof(int));
        ADDED_NODE_iSET = (int*)malloc((MAX_VERTEX_NO + 1) * 2 * sizeof(int));
        FIXED_NODE_STACK = (int*)malloc((MAX_VERTEX_NO + 1) * 2 * sizeof(int));

        iSET_State = (char*)malloc((MAX_VERTEX_NO + 1) * sizeof(char));
        iSET_Used = (char*)malloc((MAX_VERTEX_NO + 1) * sizeof(char));
        iSET_Tested = (char*)malloc((MAX_VERTEX_NO + 1) * sizeof(char));
        UNIT_STACK = (int*)malloc((MAX_VERTEX_NO + 1) * sizeof(int));
        NEW_UNIT_STACK = (int*)malloc((MAX_VERTEX_NO + 1) * sizeof(int));
        PASSIVE_iSET_STACK = (int*)malloc((MAX_VERTEX_NO + 1) * sizeof(int));
        iSET_Involved = (char*)malloc((MAX_VERTEX_NO + 1) * sizeof(char));
        CONFLICT_ISET_STACK = (int*)malloc(
            (MAX_VERTEX_NO + 1) * 10 * sizeof(int));
        REASON_STACK = (int*)malloc((MAX_VERTEX_NO + 1) * 10 * sizeof(int));
        Is_Tested = (char*)malloc((MAX_VERTEX_NO + 1) * sizeof(char));
    }
    void store_maximum_clique(int node, int print_info) {
        if (CUR_CLQ_SIZE == 0)
            PUSH_OPT(node, Clique_Stack);
        else
            PUSH_OPT(Second_Name[node], Clique_Stack);
        MAX_CLQ_SIZE = ptr(Clique_Stack);
        memcpy(MaxCLQ_Stack, Clique_Stack, MAX_CLQ_SIZE * sizeof(int));
        ptr(Candidate_Stack) = nodeNum + 1;
        ptr(Cursor_Stack) = 1;
        ptr(Clique_Stack) = 0;
        Rollback_Point = 0;
        Vertex_UB[CURSOR] = MAX_CLQ_SIZE;
        if (print_info == TRUE)
            //printf("C %4d |%7d |%8d %10d %10d %10d|%10d \n", MAX_CLQ_SIZE, CURSOR, cut_ver, cut_inc, cut_iset, cut_satz, BRANCHING_COUNT);
        total_cut_ver += cut_ver;
        cut_ver = 0;
        total_cut_inc += cut_inc;
        cut_inc = 0;
        total_cut_iset += cut_iset;
        cut_iset = 0;
        total_cut_satz += cut_satz, cut_satz = 0;
        Last_Idx = CURSOR;
        if (MAX_CLQ_SIZE == START_MAXSAT_THD)
            allocate_memory_for_maxsat();
    }
    void store_maximum_clique2(int print_info) {
        MAX_CLQ_SIZE = ptr(Clique_Stack);
        memcpy(MaxCLQ_Stack, Clique_Stack, MAX_CLQ_SIZE * sizeof(int));
        //ptr(Cursor_Stack) = 1;
        ptr(Clique_Stack) = 0;
        if (print_info == TRUE)
            //printf("C %4d |%7d |%8d %10d %10d %10d|%10d \n", MAX_CLQ_SIZE, CURSOR, cut_ver, cut_inc, cut_iset, cut_satz, BRANCHING_COUNT);
        total_cut_ver += cut_ver;
        cut_ver = 0;
        total_cut_inc += cut_inc;
        cut_inc = 0;
        total_cut_iset += cut_iset;
        cut_iset = 0;
        total_cut_satz += cut_satz, cut_satz = 0;
        Last_Idx = CURSOR;
        if (MAX_CLQ_SIZE == START_MAXSAT_THD)
            allocate_memory_for_maxsat();
    }
    int reduce_subgraph(int start) {
        int* degree_counter, * where;
        int end, p1, i, node = NONE, neibor, * neibors, t, j, h, k;
        int max_degree = 0;
        N1_0 += NB_CANDIDATE;
        L1++;
        compute_subgraph_degree(start);
        D2 += Dt;
        where = Candidate_Stack + ptr(Candidate_Stack) + 1;
        degree_counter = Vertex_UB + ptr(Candidate_Stack) + 1;
        memset(degree_counter, 0, (Max_Degree + 1) * sizeof(int));

        for (node = Candidate_Stack[i = start]; node != DELIMITER; node =
            Candidate_Stack[++i]) {
            degree_counter[Node_Degree[node]]++;
            Vertex_UB[i] = node;
        }
        Vertex_UB[i] = DELIMITER;

        end = i - 1;
        j = start;
        for (i = 0; i <= Max_Degree; i++) {
            k = degree_counter[i];
            degree_counter[i] = j;
            j += k;
        }

        for (node = Vertex_UB[i = start]; node != DELIMITER; node =
            Vertex_UB[++i]) {
            t = degree_counter[Node_Degree[node]]++;
            Candidate_Stack[t] = node;
            where[node] = t;
        }

        for (i = Max_Degree; i > 0; i--) {
            degree_counter[i] = degree_counter[i - 1];
        }
        degree_counter[0] = start;
        //return FALSE;
        p1 = start;
        while ((node = Candidate_Stack[p1]) != DELIMITER) {
            if (Node_Degree[node] > max_degree)
                max_degree = Node_Degree[node];
            if (p1 < end
                && Node_Degree[node] == Node_Degree[Candidate_Stack[p1 + 1]]) {
                degree_counter[Node_Degree[node]] = p1 + 1;
            }
            if (Node_Degree[node] == end - p1) {
                if (end - p1 + 1 == MAX_CLQ_SIZE) {
                    ptr(Clique_Stack) = 0;
                    PUSH_OPT(Candidate_Stack[CURSOR], Clique_Stack);
                    while ((node = Candidate_Stack[p1++]) != DELIMITER)
                        PUSH_OPT(node, Clique_Stack);
                    store_maximum_clique2(TRUE);
                    return TRUE;
                }
                else {
                    while ((node = Candidate_Stack[++p1]) != DELIMITER)
                        Node_Degree[node] = end - p1;
                }
                break;
            }

            neibors = iSET[iSET_Index[node]];
            for (neibor = *neibors; neibor != NONE; neibor = *(++neibors))
                if (where[neibor] > p1) {
                    t = where[neibor];
                    h = degree_counter[Node_Degree[neibor]];
                    k = Candidate_Stack[h];

                    Candidate_Stack[h] = neibor;
                    where[neibor] = h;
                    Candidate_Stack[t] = k;
                    where[k] = t;

                    degree_counter[Node_Degree[neibor]]++;
                    Node_Degree[neibor]--;
                    if (Node_Degree[neibor]
                        != Node_Degree[Candidate_Stack[h - 1]]) {
                        degree_counter[Node_Degree[neibor]] = h;
                    }
                }

            p1++;
        }
        if (max_degree + 1 >= MAX_CLQ_SIZE) {
            CUR_MAX_NODE = 0;
            for (node = Candidate_Stack[i = start]; node != DELIMITER; node =
                Candidate_Stack[++i]) {
                if (Node_Degree[node] + 1 >= MAX_CLQ_SIZE)
                    break;
            }
            j = start;
            for (node = Candidate_Stack[i]; node != DELIMITER; node =
                Candidate_Stack[++i]) {
                Vertex_UB[j] = Node_Degree[node] + 1;
                Candidate_Stack[j++] = node;
                if (node > CUR_MAX_NODE)
                    CUR_MAX_NODE = node;
                N1_1++;
            }
            Candidate_Stack[j] = DELIMITER;
            ptr(Candidate_Stack) = j + 1;
            //		if (j - start > 0) {
            //			SubNB++;
            //			//N1_0 += NB_CANDIDATE;
            //			//N1_1 += j - start;
            //			//D2 += Dt;
            //			R1 += (double) (j - start) / NB_CANDIDATE;
            //			D3 += compute_subgraph_density(start);
            //		}
            return FALSE;
        }
        else {
            Vertex_UB[CURSOR] = max_degree + 2;
            return TRUE;
        }
    }

    void rebuild_matrix(int start) {
        int i = start, j = 1, node, neibor, * neibors;
        //	for (node = 1; node <= NB_NODE; node++)
        //		Node_State[node] = PASSIVE;
        for (node = Candidate_Stack[i]; node != DELIMITER; node =
            Candidate_Stack[++i]) {
            Candidate_Stack[i] = j;
            Second_Name[j] = node;
            Node_Degree[node] = j++;
            Node_State[node] = ACTIVE;
        }
        memset(Adj_Matrix, 0,
            (MAX_VERTEX_NO + 1) * (MATRIX_ROW_WIDTH) * sizeof(char));
        i = start;
        for (node = Candidate_Stack[i]; node != DELIMITER; node =
            Candidate_Stack[++i]) {
            neibors = Node_Neibors[Second_Name[node]];
            for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                if (Node_State[neibor] == ACTIVE) {
                    SET_EDGE(node, Node_Degree[neibor]);
                    SET_EDGE(Node_Degree[neibor], node);
                }
            }
            //	Node_State[Second_Name[node]] = PASSIVE;
        }
        i = start;
        for (node = Candidate_Stack[i]; node != DELIMITER; node =
            Candidate_Stack[++i]) {
            Node_State[Second_Name[node]] = PASSIVE;
        }
    }

    int cut_by_inc_ub() {
        int i = CURSOR, neibor, max = 0, * neibors;
        int node = Candidate_Stack[CURSOR];
        int start = ptr(Candidate_Stack);
        unsigned char* adj_list;
        NB_CANDIDATE = 0;
        CUR_MAX_NODE = 0;
        if (CUR_CLQ_SIZE == 0) {
            neibors = Node_Neibors[node];
            CUR_MAX_NODE = *(neibors);
            for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                i = nodeNum - neibor;
                if (max < Vertex_UB[i])
                    max = Vertex_UB[i];
                Vertex_UB[ptr(Candidate_Stack)] = Vertex_UB[i];
                PUSH_OPT(neibor, Candidate_Stack);
                NB_CANDIDATE++;
            }
        }
        else {
            adj_list = iMatrix(node);
            while (Candidate_Stack[i] != DELIMITER)i--;
            for (neibor = Candidate_Stack[++i]; neibor != DELIMITER; neibor =
                Candidate_Stack[++i]) {
                if (neibor > 0 && Matrix(adj_list, neibor) > 0) {
                    if (max < Vertex_UB[i])
                        max = Vertex_UB[i];
                    Vertex_UB[ptr(Candidate_Stack)] = Vertex_UB[i];
                    PUSH_OPT(neibor, Candidate_Stack);
                    NB_CANDIDATE++;
                }
            }
        }
        PUSH_OPT(DELIMITER, Candidate_Stack);

        /*max=NB_CANDIDATE;//disable incremental upper bound*/
        if (NB_CANDIDATE < max) {
            max = NB_CANDIDATE;
        }
        if (Vertex_UB[CURSOR] - 1 < max) {
            max = Vertex_UB[CURSOR] - 1;
        }
        if (max < MAX_CLQ_SIZE - CUR_CLQ_SIZE) {
            Vertex_UB[CURSOR] = max + 1;
            cut_inc++;
            return TRUE;
        }
        else if (CUR_CLQ_SIZE == 0) {
            if (NB_CANDIDATE < 10 && cut_by_iset_less_vertices() == TRUE) {
                Vertex_UB[CURSOR] = iSET_COUNT + 1;
                cut_iset++;
                return TRUE;
            }
            else if (reduce_subgraph(start)) {
                cut_iset++;
                return TRUE;
            }
            else if (REBUILD_MATRIX == TRUE || CUR_MAX_NODE > MAX_VERTEX_NO) {
                rebuild_matrix(start);
                REBUILD_MATRIX = TRUE;
                return FALSE;
            }
            return FALSE;
        }
        else {
            return FALSE;
        }
    }

    int find_3_clique(int node) {
        int neibor1, neibor2, neibor3, * neibors1, * neibors2, * neibors3;
        neibors1 = Node_Neibors[node];
        for (neibor1 = *neibors1; neibor1 != NONE; neibor1 = *(++neibors1)) {
            neibors2 = neibors1 + 1;
            for (neibor2 = *neibors2; neibor2 != NONE; neibor2 = *(++neibors2)) {
                if (is_adjacent(neibor1, neibor2) == TRUE) {
                    neibors3 = neibors2 + 1;
                    for (neibor3 = *neibors3; neibor3 != NONE; neibor3 =
                        *(++neibors3)) {
                        if (is_adjacent(neibor1, neibor3) == TRUE
                            && is_adjacent(neibor2, neibor3) == TRUE) {
                            MaxCLQ_Stack[0] = Old_Name[node];
                            MaxCLQ_Stack[1] = Old_Name[neibor1];
                            MaxCLQ_Stack[2] = Old_Name[neibor2];
                            MaxCLQ_Stack[3] = Old_Name[neibor3];
                            MAX_CLQ_SIZE = 4;
                            return TRUE;
                        }
                    }
                }
            }
        }
        return FALSE;
    }

    void init_for_search(int using_init_clique) {
        int i, node;
        int neibor, neibor2, * neibors, * neibors2;
        cut_ver = 0;
        cut_inc = 0;
        cut_iset = 0;
        cut_satz = 0;
        total_cut_ver = 0;
        total_cut_inc = 0;
        total_cut_iset = 0;
        total_cut_satz = 0;

        Branches_Nodes[0] = 0;
        Branches_Nodes[1] = 0;
        Branches_Nodes[2] = 0;
        Branches_Nodes[3] = 0;
        Branches_Nodes[4] = 0;
        Branches_Nodes[5] = 0;

        Last_Idx = nodeNum;
        //NB_BACK_CLIQUE = 0;
        MAX_CLQ_SIZE = 0;
        ptr(Clique_Stack) = 0;
        ptr(Cursor_Stack) = 0;
        Rollback_Point = 0;
        PUSH_OPT(nodeNum - INIT_CLQ_SIZE - 1, Cursor_Stack);
        MAX_CLQ_SIZE = INIT_CLQ_SIZE;
        for (i = 0; i < ptr(Candidate_Stack) - 1; i++) {
            node = Candidate_Stack[i];
            Vertex_UB[i] = Node_Degree[node] + 1;
            if (INIT_CLQ_SIZE == 3 && Vertex_UB[i] > 3) {
                if (find_3_clique(node) == TRUE) {
                    //printf("I improve the initial clique to 4\n");
                    INIT_CLQ_SIZE = 4;
                }
                else {
                    Vertex_UB[i] = 3;
                }
            }
            if (INIT_CLQ_SIZE == 4 && Vertex_UB[i] == 5) {
                neibors = Node_Neibors[node];
                for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                    neibors2 = neibors + 1;
                    for (neibor2 = *neibors2; neibor2 != NONE; neibor2 =
                        *(++neibors2)) {
                        if (is_adjacent(neibor, neibor2) == FALSE) {
                            Vertex_UB[i] = 4;
                            break;
                        }
                    }
                    if (Vertex_UB[i] == 4)
                        break;
                }
            }
        }
    }

    void allocate_memory() {
        int i;
        Second_Name = (int*)malloc((MAX_VERTEX_NO + 1) * sizeof(int));
        iSET = (int**)malloc((MAX_VERTEX_NO + 1) * sizeof(int*));
        iSET[0] = (int*)malloc(
            (MAX_VERTEX_NO + 1) * (MAX_VERTEX_NO + 1) * sizeof(int));
        for (i = 1; i < MAX_VERTEX_NO; i++) {
            iSET[i] = iSET[i - 1] + MAX_VERTEX_NO + 1;
        }
        iSET_Size = (int*)malloc((MAX_VERTEX_NO + 1) * sizeof(int));
        iSET_Index = (int*)malloc((nodeNum + 1) * sizeof(int));

        if (INIT_CLQ_SIZE >= START_MAXSAT_THD)
            allocate_memory_for_maxsat();
    }

    bool search_maxclique(int cutoff, int using_init_clique, double MaxTime) {
        int node;
        init_for_search(using_init_clique);
        BRANCHING_COUNT = 0;
        if (using_init_clique == TRUE) {
            //printf("C  -----------------------------------------------------------------\n");
            //printf( "C  Size|   Index|NB_Vertex  NB_IncUB    NB_Iset  NB_MaxSat|  NB_Branch\n");
        }
        clock_t start = clock();
        while (CURSOR > 0 && (clock() - start) * 1.0 / CLOCKS_PER_SEC < MaxTime) {
            node = Candidate_Stack[--CURSOR];
            if (CUR_CLQ_SIZE > 0 && node > 0)
                continue;
            if (node == DELIMITER) {
                ptr(Candidate_Stack) = CURSOR + 1;
                ptr(Cursor_Stack)--;
                ptr(Clique_Stack)--;
                Vertex_UB[CURSOR] = MAX_CLQ_SIZE - CUR_CLQ_SIZE;
            }
            else {
                if (node < 0) {
                    node = -node;
                    Candidate_Stack[CURSOR] = -Candidate_Stack[CURSOR];
                }
                if (MAX_CLQ_SIZE == CUR_CLQ_SIZE) {
                    store_maximum_clique(node, using_init_clique);
                }
                else if (Vertex_UB[CURSOR] <= MAX_CLQ_SIZE - CUR_CLQ_SIZE) {
                    cut_ver++;
                }
                else {
                    BRANCHING_COUNT++;
                    Rollback_Point = ptr(Candidate_Stack);

                    if (cut_by_inc_ub() == TRUE

                        || cut_by_iset_last_renumber() == TRUE

                        || (MAX_CLQ_SIZE >= START_MAXSAT_THD && cut_by_inc_maxsat_eliminate_first() == TRUE)
                        ) {
                        ptr(Candidate_Stack) = Rollback_Point;
                    }
                    else {
                        if (CUR_CLQ_SIZE == 0)
                            PUSH_OPT(node, Clique_Stack);
                        else
                            PUSH_OPT(Second_Name[node], Clique_Stack);
                        PUSH_OPT(Branching_Point, Cursor_Stack);
                    }
                }
            }
        }
        if (using_init_clique == TRUE) {
            //printf("C  ----------------------------------------------------------------\n");
            //printf("C %4d |%7d |%8d %10d %10d %10d|%10d \n", MAX_CLQ_SIZE, CURSOR, cut_ver, cut_inc, cut_iset, cut_satz, BRANCHING_COUNT);
            total_cut_ver += cut_ver;
            total_cut_inc += cut_inc;
            total_cut_iset += cut_iset;
            total_cut_satz += cut_satz;
            //printf("C  ----------------------------------------------------------------\n");
            //printf("C %4d |%7d |%8d %10d %10d %10d|%10d \n", MAX_CLQ_SIZE, CURSOR, total_cut_ver, total_cut_inc, total_cut_iset, total_cut_satz, BRANCHING_COUNT);
        }
        if (CURSOR > 0) 
        { 
            return false; 
        }
        else return true;
        
    }

    void printMaxClique(Vec<int>& solution) {
        int i;
        //printf("M ");
        if (INIT_CLQ_SIZE < MAX_CLQ_SIZE) {
            for (i = 0; i < MAX_CLQ_SIZE; i++) {
                //printf("%d ", Old_Name[MaxCLQ_Stack[i]]);
                solution.push_back(Old_Name[MaxCLQ_Stack[i]]);
            }

        }
        else {
            for (i = 0; i < MAX_CLQ_SIZE; i++) {
                //printf("%d ", MaxCLQ_Stack[i]);
                solution.push_back(MaxCLQ_Stack[i]);
            }

        }
        //printf("\n");
    }
    void build_init_matrix() {
        int node, neibor, * neibors;
        MATRIX_ROW_WIDTH = MAX_VERTEX_NO / 8 + 1;
        Adj_Matrix = (unsigned char*)malloc(
            (MAX_VERTEX_NO + 1) * MATRIX_ROW_WIDTH * sizeof(char));

        memset(Adj_Matrix, 0,
            (MAX_VERTEX_NO + 1) * MATRIX_ROW_WIDTH * sizeof(char));

        for (node = 1; node <= MAX_VERTEX_NO; node++) {
            Second_Name[node] = node;
            neibors = Node_Neibors[node];
            for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                SET_EDGE(node, neibor);
                SET_EDGE(neibor, node);
            }
        }
    }
    int* Adj_List;
#define New_Name Node_Degree
    int search_in_2_k_core_graph() {
        int i = 0, node, neibor1, neibor2, neibor, * neibors, find = FALSE;
        for (i = 0; i < nodeNum; i++) {
            node = Candidate_Stack[i];
            if (Node_Degree[node] == 2) {
                neibor1 = Node_Neibors[node][0];
                neibor2 = Node_Neibors[node][1];

                neibors = Node_Neibors[neibor1];
                for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                    if (neibor == neibor2) {
                        find = TRUE;
                        break;
                    }
                }
                if (find == TRUE) {
                    MaxCLQ_Stack[0] = Old_Name[node];
                    MaxCLQ_Stack[1] = Old_Name[neibor1];
                    MaxCLQ_Stack[2] = Old_Name[neibor2];
                    INIT_CLQ_SIZE = 3;
                    MAX_CLQ_SIZE = 3;
                    break;
                }
            }
            else if (Node_Degree[node] > 2) {
                break;
            }
        }
        if (find == TRUE) {
#if LOCALTEST
            //printf("I find maximum clique in initial phase.");
#endif
            return TRUE;
        }
        else if (i == nodeNum) {
#if LOCALTEST
            //printf("I find maximum clique in initial phase.");
#endif
            MAX_CLQ_SIZE = 2;
            return TRUE;
        }
        else {
#if LOCALTEST 
            //printf("I not find maximum clique in initial phase!!");
#endif
            return FALSE;
        }
    }

    void free_block() {
        int i = 0;
        for (i = 0; i < BLOCK_COUNT; i++)
            free(BLOCK_LIST[i]);
    }

    void reduce_instance() {
        int i, j, nb, node, * neibors, * neibors2, * addr;
        MAX_VERTEX_NO = 0;
        j = 0;
        for (i = 0; i < nodeNum; i++) {
            node = Candidate_Stack[i];
            if (CORE_NO[i] < INIT_CLQ_SIZE) {
                Node_State[node] = PASSIVE;
            }
            else {
                Candidate_Stack[j++] = node;
                Node_State[node] = ACTIVE;
            }
        }
        nodeNum = j;
        N0_1 = nodeNum;
        Candidate_Stack[j] = DELIMITER;
        ptr(Candidate_Stack) = j + 1;
        Old_Name = (int*)malloc((nodeNum + 1) * sizeof(int));
        for (i = 0; i < nodeNum; i++) {
            Old_Name[nodeNum - i] = Candidate_Stack[i];
            New_Name[Candidate_Stack[i]] = nodeNum - i;
            Candidate_Stack[i] = nodeNum - i;
        }

        edgeNum = 0;
        for (i = nodeNum; i > 0; i--) {
            neibors = Node_Neibors[Old_Name[i]];
            neibors2 = neibors;
            nb = 0;
            for (node = *neibors; node != NONE; node = *(++neibors)) {
                if (Node_State[node] == ACTIVE && New_Name[node] < i) {
                    (*neibors2) = New_Name[node];
                    neibors2++;
                    nb++;
                }
            }
            (*neibors2) = NONE;
            edgeNum += nb;
            qsort(Node_Neibors[Old_Name[i]], nb, sizeof(int), int_cmp);
        }

        Adj_List = (int*)malloc((edgeNum + nodeNum) * sizeof(int));
        addr = Adj_List;

        if (Adj_List == NULL) {
            //printf("can't allocate memory for Adj_List!\n");
            exit(0);
        }

        for (i = nodeNum; i > 0; i--) {
            Node_Degree[i] = 0;
            Node_State[i] = PASSIVE;
            neibors = Node_Neibors[Old_Name[i]];
            for (node = *neibors; node != NONE; node = *(++neibors)) {
                *(addr++) = node;
                Node_Degree[i]++;
            }
            *(addr++) = NONE;
            if (Node_Degree[i] > MAX_VERTEX_NO)
                MAX_VERTEX_NO = Node_Degree[i];
        }
        //free(Init_Adj_List);
        free_block();
        Node_Neibors[nodeNum] = Adj_List;
        for (i = nodeNum - 1; i > 0; i--) {
            Node_Neibors[i] = Node_Neibors[i + 1] + Node_Degree[i + 1] + 1;
        }
        D1 = ((float)edgeNum * 2 / nodeNum / (nodeNum - 1));
        //printf("I the reduced graph #node %d #edge %d #density %9.8f\n", nodeNum,edgeNum, ((float)edgeNum * 2 / nodeNum / (nodeNum - 1)));
        //printf("I the largest subgraph is %d\n", MAX_VERTEX_NO);

    }

    int initialize(int ordering) {
        clock_t state = clock();
        int r = sort_by_degeneracy_ordering();
        if (r == FALSE) {
            //printf("*****************\n");
            reduce_instance();
            
            if (K_CORE_G <= 2) {
                r = search_in_2_k_core_graph();
             
            }
        }
        INIT_TIME = (clock() - state) * 1.0 / CLOCKS_PER_SEC;

        //printf("I the initial time is %4.2lf \n", INIT_TIME);
       //if(r==0) //printf("init!");
        return r;
    }

    void print_branches() {
        int i = 0;
        for (i = 0; i <= nodeNum; i++) {
            //printf("%3d", Branches[i]);
        }
        //printf("\n");
    }
    void print_version() {
        //printf("# Hello! I am LMC with Incremental MaxSAT Reasoning and 2 Level Cores Decomposition. build at %s %s.\n", __TIME__, __DATE__);
        return;
    }

    char* getInstanceName(char* s) {
        if (strrchr(s, '/') == NULL)
            return s;
        else
            return strrchr(s, '/') + 1;
    }
    bool InstanceLoad(Vec<NewEdge>& newEdges, int newNodeNum) {
        int max_node = NONE, offset = 0;
        int nb_edge = 0;
        memset(Node_Degree, 0, newNodeNum * sizeof(int));
        for (auto i = newEdges.begin(); i != newEdges.end(); ++i) {
            Node_Degree[(*i).src]++;
            Node_Degree[(*i).dst]++;
            if ((*i).dst > max_node)
                max_node = (*i).dst;
            if ((*i).src > max_node)
                max_node = (*i).src;
            if ((*i).src == 0 || (*i).dst == 0) {
                offset = 1;
            }
            nb_edge++;
        }
        nodeNum = max_node;
        nodeNum = nodeNum + offset;

        Node_Neibors = (int**)malloc((nodeNum + 1) * sizeof(int*));
        allcoate_memory_for_adjacency_list(nodeNum, nb_edge, offset);
        memset(Node_Degree, 0, (nodeNum + 1) * sizeof(int));
        nb_edge = 0;
        for (auto i = newEdges.begin(); i != newEdges.end(); ++i) {
            int l_node = (*i).src, r_node = (*i).dst;
            if (l_node >= 0 && r_node >= 0 && l_node != r_node) {
                l_node += offset;
                r_node += offset;
                int j = 0;
                for (j = 0; j < Node_Degree[l_node]; j++) {
                    if (Node_Neibors[l_node][j] == r_node)
                        break;
                }
                if (j == Node_Degree[l_node]) {
                    Node_Neibors[l_node][Node_Degree[l_node]] = r_node;
                    Node_Neibors[r_node][Node_Degree[r_node]] = l_node;
                    Node_Degree[l_node]++;
                    Node_Degree[r_node]++;
                    nb_edge++;
                }
            }
        }
        edgeNum = nb_edge;
        Max_Degree = 0;
        N0_0 = nodeNum;
        for (int node = 1; node <= nodeNum; node++) {
            Node_Neibors[node][Node_Degree[node]] = NONE;
            if (Node_Degree[node] > Max_Degree)
                Max_Degree = Node_Degree[node];
        }
        NB_NODE_O = nodeNum;
        NB_EDGE_O = edgeNum;
        D0 = ((float)edgeNum * 2 / nodeNum / (nodeNum - 1));
        return true;
    }
    bool loadOriginresult(Vec<NewEdge>& newEdges, int format) {
        int j, l_node, r_node, nb_edge = 0, max_node = NONE, offset = 0;
        int node = 1;
        char line[128], word[10];
        memset(Node_Degree, 0, MAX_NODE * sizeof(int));
        for (auto i = newEdges.begin(); i != newEdges.end(); ++i) {
            Node_Degree[(*i).src]++;
            Node_Degree[(*i).dst]++;
            if ((*i).dst > max_node)
                max_node = (*i).dst;
            if ((*i).src == 0 || (*i).dst == 0) {
                offset = 1;
            }
        }
        nodeNum = max_node;
        nodeNum = nodeNum + offset;
        if (nodeNum > MAX_NODE) {
            exit(0);
        }

        Node_Neibors = (int**)malloc((nodeNum + 1) * sizeof(int*));
        allcoate_memory_for_adjacency_list(nodeNum, nb_edge, offset);
        memset(Node_Degree, 0, (nodeNum + 1) * sizeof(int));

        nb_edge = 0;
        for (auto i = newEdges.begin(); i != newEdges.end(); ++i) {
            int l_node = (*i).src, r_node = (*i).dst;
            if (l_node >= 0 && r_node >= 0 && l_node != r_node) {
                l_node += offset;
                r_node += offset;
                int j = 0;
                for (j = 0; j < Node_Degree[l_node]; j++) {
                    if (Node_Neibors[l_node][j] == r_node)
                        break;
                }
                if (j == Node_Degree[l_node]) {
                    Node_Neibors[l_node][Node_Degree[l_node]] = r_node;
                    Node_Neibors[r_node][Node_Degree[r_node]] = l_node;
                    Node_Degree[l_node]++;
                    Node_Degree[r_node]++;
                    nb_edge++;
                }
            }
        }
        edgeNum = nb_edge;
        Max_Degree = 0;
        N0_0 = nodeNum;
        for (node = 1; node <= nodeNum; node++) {
            Node_Neibors[node][Node_Degree[node]] = NONE;

            if (Node_Degree[node] > Max_Degree)
                Max_Degree = Node_Degree[node];
        }
        return true;
    }
    void check_result(Vec<NewEdge>& newEdges, int* solvers) {

        if (FORMAT == 1) {
            loadOriginresult(newEdges, 1);
        }
        else if (FORMAT == 2) {
            loadOriginresult(newEdges, 2);
        }
        else {
            loadOriginresult(newEdges, 1);
        }
        /*//printf("R Instance Information: #node=%d #edge=%d density=%9.8f\n", nodeNum,
            edgeNum, ((float)edgeNum * 2 / nodeNum / (nodeNum - 1)));*/
        check_clique_in_result_file(solvers);

    }
    void minVertexCovering(Vec<bool>& outputSolver, Vec<int>& maxCliqueSolution) {
       // printf("%d:\n", MAX_CLQ_SIZE);
        if (INIT_CLQ_SIZE < MAX_CLQ_SIZE) {
            for (int i = 0; i < MAX_CLQ_SIZE; i++) {
                outputSolver[Old_Name[MaxCLQ_Stack[i]]] = 1;
                //printf("%d ", Old_Name[MaxCLQ_Stack[i]]);
                maxCliqueSolution.push_back(Old_Name[MaxCLQ_Stack[i]]);
            }
        }
        else {
            for (int i = 0; i < MAX_CLQ_SIZE; i++) {
                outputSolver[MaxCLQ_Stack[i]] = 1;
                //printf("%d ", MaxCLQ_Stack[i]);
                maxCliqueSolution.push_back(MaxCLQ_Stack[i]);
            }
        }
    }
    bool Checker(Vec<int>& solution, Vec<NewEdge>& newEdges, int newNodeNum) {
        InstanceLoad(newEdges, newNodeNum);
        int clique_size = solution.size();
        int i = 0;
        if (clique_size > 0) {
            for (i = 0; i < clique_size; i++) {
                for (int j = i + 1; j < clique_size; j++) {
                    if (is_adjacent(solution[i], solution[j]) == false) {
                        return false;
                    }
                }
            }
            ////printf("#true\n");
        }
        else {
            ////printf("#NONE\n");
        }
        return true;
    }
public:
    bool maxCliqueSolver(Vec<NewEdge>& newEdges, int newNodeNum, Vec<bool>& Solver, double maxTime) {
        int ordering = -1;
        double MaxTime = maxTime;
        clock_t state = clock();
        bool flag = true;
        Vec<int> maxCliqueSolution;
        int nodeNum_O, edgeNum_O;
        if (InstanceLoad(newEdges, newNodeNum)) {
            nodeNum_O = nodeNum;
            edgeNum_O = edgeNum;

            if (initialize(ordering) == false) {
                allocate_memory();
                build_init_matrix();
                flag = search_maxclique(0, true, MaxTime);
            }
            
            if (flag == false) {
                return false;
            }
            minVertexCovering(Solver, maxCliqueSolution);
            if (!Checker(maxCliqueSolution, newEdges, newNodeNum)) {
                return false;
           }
        }
        return true;
    }
};