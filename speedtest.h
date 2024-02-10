#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>

#define WARMUP_QUERIES 100000000

typedef uint64_t (*FUNC_PTR)(uint64_t);

// mod is num_ones for select and num_bits for rank
// rank_or_select = 1 if select and 0 if rank
int fp_speed_test(FUNC_PTR fun, uint64_t num_queries, uint64_t mod, uint64_t rank_or_select) {
    
    srand(time(NULL));

    clock_t start, end;
    double avg_time = 0;
    if (num_queries < WARMUP_QUERIES) {
        num_queries = WARMUP_QUERIES;
    }

    uint64_t *numbers = (uint64_t *)malloc(sizeof(uint64_t) * num_queries);
    uint64_t *numbers1 = (uint64_t *)malloc(sizeof(uint64_t) * num_queries);

    for (int i = 0; i < num_queries; i++) {
        uint64_t my_rand = rand();
        my_rand = ((my_rand << 32) | rand()) % (mod);
        numbers[i] = my_rand + rank_or_select;
    }

    for (int i = 0; i < num_queries; i++) {
        uint64_t my_rand = rand();
        my_rand = ((my_rand << 32) | rand()) % (mod);
        numbers1[i] = my_rand + rank_or_select;
    }

    volatile uint64_t result;
    // WARMUP
    for (int i = 0; i < WARMUP_QUERIES; i++)
        result = fun(numbers[i]);

    free(numbers);
    // TEST
    start = clock();
    for (int i = 0; i < num_queries; i++)
        result = fun(numbers1[i]);
    

    end = clock();
    clock_t clock_time = end - start;
    // takes average and converts to nano seconds
    avg_time = (((((double)clock_time) * 1000000000)) / num_queries) / CLOCKS_PER_SEC;
    printf("%lf\n", avg_time);
    free(numbers1);

    return 0;
}