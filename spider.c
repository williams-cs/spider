// use aligned malloc
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <x86intrin.h>
#include "helper.h"

// SECTION - Basic Block Queries
static inline uint64_t rank_512(uint64_t *bits, uint64_t curr_uint, uint64_t nbits) {
    uint64_t last_uint = (nbits) / 64llu;

    uint64_t pop_val = 0;
    uint64_t final;
    if ((last_uint)) {

        pop_val += __builtin_popcountll(bits[curr_uint] & BIT_MASK);

        for (int i = 1; i < last_uint; i++) {
            pop_val += __builtin_popcountll(bits[curr_uint + i]);
        }

        final = bits[curr_uint + last_uint] >> (63 - ((nbits) & 63llu));
    }
    else {
        final = (bits[curr_uint] & BIT_MASK) >> (63 - ((nbits)));
    }
    pop_val += __builtin_popcountll(final);
    return pop_val;
}

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
//!SECTION

//SECTION - Quieries

//SECTION - Rank
static inline uint64_t rank(uint64_t pos) {
    uint64_t section = pos / 496;
    uint64_t l0_i = section >> 7;
    uint64_t start = section << 3;

    return l0_a[l0_i] + ((modified_bit_array[start] & MARK_MASK) >> 48) + rank_512(modified_bit_array, start, (pos - (section * 496)) + 16);
}
//!SECTION

//SECTION - Select
static inline uint64_t select_hl(uint64_t ones) {

    uint64_t l0_guess = hl_select_array[((ones -1) >> log_hl_select_width)];

    while (ones <= l0_a[l0_guess]) {

        l0_guess--;
    }

    while (ones > l0_a[l0_guess + 1]) {

        l0_guess++;
    }

    uint64_t select_i = (ones - 1) >> log_ones_per_slot;
    uint64_t s_curr = ll_select_array[select_i]; 
    uint64_t s_next = ll_select_array[select_i + 1];

    // handles the case if our select array spans 2 superblocks
    // must handle if the target is in the first superblock and if it is second
    if (s_next < s_curr) {
        if ((l0_a[l0_guess] >> log_ones_per_slot) == select_i) {

            s_curr -= SUPERBLOCK_SIZE;
        } else {

            s_next += SUPERBLOCK_SIZE;
        }
    }

    uint64_t curr512_i = ((((ones - (select_i * ones_per_slot)) * (s_next - s_curr) >> log_ones_per_slot)) + (s_curr + l0_guess * SUPERBLOCK_SIZE)) / 496;
    

    while (ones <= ((modified_bit_array[curr512_i * 8] & MARK_MASK) >> 48) + l0_a[curr512_i >> 7]) {

        curr512_i--;
    }

    while (ones > ((((modified_bit_array[(curr512_i + 1) * 8] & MARK_MASK)) >> 48) + l0_a[(curr512_i+1) >> 7])) {

        curr512_i++;
    }

    return (curr512_i * 496) + select_512(curr512_i * 8, ones - (((modified_bit_array[curr512_i * 8] & MARK_MASK) >> 48) + l0_a[curr512_i >> 7]));
}
//!SECTION


//SECTION - Driver
int main(int argc, char *argv[]) {
    uint64_t query;

    if (argc < 4) {
        help();
    }

    clock_t b_start = clock();
    bit_meta pack = modify_package_byte_file(argv[1], atoll(argv[2]), 4);
    uint64_t total_size = (pack.num_elements * 64) - pack.nbits + (pack.l0_size * 64);
    correctness_meta tester = {};


    if (argv[3][0] == 's' || argv[3][0] == 'd' || argv[3][0] == 'q') {

        total_size += build_select_from_modified(pack);
    }
    clock_t b_end = clock();
    #ifndef VERBOSE
        printf("%lf\n", (double)(b_end - b_start) / (double)CLOCKS_PER_SEC);
        printf("Percent Total Size: %lf\n", ((double)total_size * 100) / (double)(pack.nbits));
    #endif
    

    switch(argv[3][0]) {
        case 'r':
            fp_speed_test(&rank, WARMUP_QUERIES, pack.nbits);
            break;
        case 's':
            fp_speed_test(&select_hl, WARMUP_QUERIES, pack.num_ones);
            break;
        
        case 'c':
            printf("Running Rank Correctness:\n");
            tester = package_byte_file_for_correctness(argv[1], atoll(argv[2]));
            fp_correct_test(&rank, tester, 100000000000lu, 1);
            break;

        case 'd':
            printf("Running Select Correctness:\n");
            tester = package_byte_file_for_correctness(argv[1], atoll(argv[2]));
            fp_correct_test(&select_hl, tester, 100000000000lu, 0);
            break;

        case 'q':

            if (argv[3][1] == 'r') {
                printf("Rank Query on Range: 0-%ld\n", pack.nbits - 1);
                printf("Enter Query or -1 to quit: ");
                if (scanf("%ld", &query) < 1) {
                    printf("Must provide a number\n");
                }
                while (query != -1) {
                    if (query < 0 || query >= pack.nbits) {
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
        case '?':
            help();
            break;
    }

    free_metadata();

}
//!SECTION
