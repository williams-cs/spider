#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <x86intrin.h>

#include "helper.h"
#include "speedtest.h"

static inline uint64_t select_512(uint64_t starting_offset, uint64_t ones) {
   int curr_uint = 0;
   uint64_t to_pop64 = ((BIT_MASK & modified_bit_array[(starting_offset)]) << 16);
   int pop64 = __builtin_popcountll(to_pop64);
   int adjust = 0; 
   

   while (ones > pop64 && curr_uint < 7) {
        adjust = 16;
         ones -= pop64;
         curr_uint++;
         to_pop64 = (modified_bit_array[(starting_offset) + curr_uint]);
         pop64 = __builtin_popcountll(to_pop64);
   }

   return curr_uint * 64 + _lzcnt_u64(_pdep_u64(1lu << (__builtin_popcountll(to_pop64) - ones), to_pop64)) - adjust;
}

uint64_t select_hl(uint64_t ones) {

    uint64_t l0_guess = select_array_1L[ones >> log_ones_per_slot] / SUPERBLOCK_SIZE;


    uint64_t select_i = (ones - 1) >> log_ones_per_slot;

    uint64_t s_curr = select_array_1L[select_i]; 
    uint64_t s_next = select_array_1L[select_i + 1];

    uint64_t curr512_i = (((((ones - ((select_i * ones_per_slot))) * (s_next - s_curr)) >> log_ones_per_slot)) + (s_curr)) / 496;


    while (ones <= ((modified_bit_array[curr512_i * 8] & MARK_MASK) >> 48) + l0_a[curr512_i >> 7]) {

        curr512_i--;
    }

    while (ones > ((((modified_bit_array[(curr512_i + 1) * 8] & MARK_MASK)) >> 48) + l0_a[(curr512_i+1) >> 7])) {
        curr512_i++;
    }

    return (curr512_i * 496) + select_512(curr512_i * 8, ones - (((modified_bit_array[curr512_i * 8] & MARK_MASK) >> 48) + l0_a[curr512_i >> 7]));

    
}




//SECTION - Main
int main(int argc, char *argv[]) {
    uint64_t query;
    clock_t b_start = clock();
    bit_meta pack = read_and_build_rank(argv[1], atoll(argv[2]), 6);
    uint64_t total_size = (pack.num_elements * 64) - pack.nbits + (pack.l0_size * 64);


    if (argv[3][0] == 's' || strcmp(argv[3], "qs") == 0) {


        total_size += build_1L_select_from_modified(pack);
    }
    clock_t b_end = clock();
    #ifndef VERBOSE
        printf("%lf\n", (double)(b_end - b_start) / (double)CLOCKS_PER_SEC);
        printf("Percent Total Size: %lf\n", ((double)total_size * 100) / (double)(pack.nbits));
    #endif
    

    switch(argv[3][0]) {
        case 'r':
            printf("Rank not supported: use ./spider instead\n");

            exit(1);
            break;
        case 's':
            fp_speed_test(&select_hl, WARMUP_QUERIES, pack.num_ones, 1);
            break;

        case 'q':

            if (argv[3][1] == 'r') {
                printf("Rank not supported: use ./spider instead\n");
                exit(1);

            }
        
            if (argv[3][1] == 's') {
                printf("Select Query on Range: 1-%ld\n", pack.num_ones);
                printf("Enter Query or -1 to quit: ");
                if (scanf("%ld", &query) < 1) {
                    printf("Must provide a number\n");
                }
                while (query != -1) {
                    if (query < 1 || query > pack.num_ones) {
                        printf("Query out of Range\n");
                    } else {
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

    free(modified_bit_array);
    free(select_array_1L);
    free(l0_a);

}
//!SECTION
