// use aligned malloc
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <x86intrin.h>
#include "ni-helper.h"
#include "speedtest.h"


static inline uint64_t rank_512(uint64_t *bits, uint64_t curr_uint, uint64_t nbits) {
    uint64_t last_uint = (nbits) / 64llu;

    uint64_t pop_val = 0;


    for (int i = 0; i < last_uint; i++) {
        pop_val += __builtin_popcountll(bits[curr_uint + i]);
    }

    uint64_t final = bits[curr_uint + last_uint] >> (63 - ((nbits) & 63llu));

    pop_val += __builtin_popcountll(final);
    return pop_val;
}

static inline uint64_t select_512(uint64_t starting_offset, uint64_t ones) {
   int curr_uint = 0;
   uint64_t to_pop64 = bit_a[(starting_offset)];
   int pop64 = __builtin_popcountll(to_pop64);
   while (ones > pop64 && curr_uint < 7) {
         ones -= pop64;
         curr_uint++;
         to_pop64 = (bit_a[(starting_offset) + curr_uint]);
         pop64 = __builtin_popcountll(to_pop64);
   }

   return curr_uint * 64 + _lzcnt_u64(_pdep_u64(1lu << (__builtin_popcountll(to_pop64) - ones), to_pop64));
}

//SECTION - Quieries

uint64_t rank(uint64_t pos) {
    uint64_t l0_i = pos >> 16;
    uint64_t l1_i = pos >> 9;
    uint64_t start = l1_i << 3;
    uint64_t so_far = l0_a[l0_i] + l1_a[l1_i];
    return so_far + rank_512(bit_a, start, pos - (l1_i * BASICBLOCK_SIZE));
}

uint64_t select_hl(uint64_t ones) {

    uint64_t l0_guess = select_a[ones >> log_ones_per_slot] / SUPERBLOCK_SIZE;

    uint64_t select_i = (ones - 1) >> log_ones_per_slot;

    uint64_t s_curr = select_a[select_i]; 
    uint64_t s_next = select_a[select_i + 1];

    uint64_t curr512_i = ((((ones - ((select_i * ones_per_slot))) * (s_next - s_curr)) >> log_ones_per_slot) + s_curr) >> 9;
    

    while (ones <= ((l1_a[curr512_i]) + l0_a[curr512_i >> 7])) {

        curr512_i--;
    }

    while (ones > (l1_a[(curr512_i + 1)] + l0_a[(curr512_i + 1) >> 7])) {

        curr512_i++;
    }

    return (curr512_i * 512) + select_512(curr512_i * 8, ones - ((l1_a[curr512_i]) + l0_a[curr512_i >> 7]));
}

//!SECTION

//SECTION - Main
int main(int argc, char *argv[]) {
    uint64_t query;
    clock_t b_start = clock();
    meta64 bit_data = safe_package_byte_file(argv[1], atoll(argv[2]));
    bit_a = bit_data.bit_array;
    uint64_t total_size = build_ni_rank(bit_data);
    
    if (argv[3][0] == 's' || strcmp(argv[3], "qs") == 0) {

        total_size += build_ni_select(bit_data);
    }

    clock_t b_end = clock();
    #ifndef VERBOSE
        printf("%lf\n", (double)(b_end - b_start) / (double)CLOCKS_PER_SEC);
        printf("Percent size: %lf\n", ((double)total_size * 100) / (double)(bit_data.nbits));
    #endif
    
    if (argc != 4) {
        help();
    }

    switch(argv[3][0]) {
        case 'r':
            fp_speed_test(&rank, WARMUP_QUERIES, bit_data.nbits, 0);
            break;
        case 's':
            fp_speed_test(&select_hl, WARMUP_QUERIES, bit_data.num_ones, 1);
            break;
            
        case 'q':

            if (argv[3][1] == 'r') {
                printf("Rank Query on Range: 0-%ld\n", bit_data.nbits - 1);
                printf("Enter Query or -1 to quit: ");
                if (scanf("%ld", &query) < 1) {
                    printf("Must provide a number\n");
                }
                while (query != -1) {
                    if (query < 0 || query >= bit_data.nbits) {
                        printf("Query out of Range\n");
                    } else {
                        printf("Rank(%ld) = %ld\n", query, rank(query));

                    }
                    printf("Enter Query or -1 to quit: ");
                    if (scanf("%ld", &query) < 1) {
                        printf("\nERROR: Must provide a number\n");
                        exit(-1);

                    }
                }

            }
        
            if (argv[3][1] == 's') {
                build_ni_select(bit_data);
                printf("Select Query on Range: 1-%ld\n", bit_data.num_ones);
                printf("Enter Query or -1 to quit: ");
                if (scanf("%ld", &query) < 1) {
                    printf("Must provide a number\n");
                }
                while (query != -1) {
                    if (query < 1 || query > bit_data.num_ones) {
                        printf("Query out of Range\n");
                    } else {
                        printf("%ld\n", query);
                        printf("Select(%ld) = %ld\n", query, select_hl(query));

                    }
                    printf("Enter Query or -1 to quit: ");
                    if (scanf("%ld", &query) < 1) {
                        printf("\nERROR: Must provide a number\n");
                        exit(-1);

                    }
                }
            }

            break;
        default:
            help();
            break;
    }

    free(bit_a);
    free(select_a);
    free(l0_a);
    free(l1_a);


}
//!SECTION
