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


#define SUPERBLOCK_SIZE 65536
#define BASICBLOCK_SIZE 512

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

uint64_t select_hl(uint64_t ones) {

    uint64_t l0_guess = (hl_select_a[(ones >> log_hl_select_width)] + 32768) / SUPERBLOCK_SIZE;

    while (ones <= l0_a[l0_guess]) {

        l0_guess--;
    }

    while (ones > l0_a[l0_guess + 1]) {

        l0_guess++;
    }

    uint64_t select_i = (ones - 1) >> log_ones_per_slot;

    uint64_t s_curr = ll_select_a[select_i]; 
    uint64_t s_next = ll_select_a[select_i + 1];

    // handles the case if our select array spans 2 superblocks
    // must handle if the target is in the first superblock and if it is second
    if (s_next < s_curr) {
        
        if ((l0_a[l0_guess] >> log_ones_per_slot) == select_i) {
            s_curr -= SUPERBLOCK_SIZE;

        } else {
            s_next += SUPERBLOCK_SIZE;
        }
    }

    uint64_t curr512_i = ((((ones - ((select_i * ones_per_slot))) * (s_next - s_curr)) >> log_ones_per_slot) + (((s_curr + l0_guess * SUPERBLOCK_SIZE)))) >> 9;


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
    double float_ratio = ((double)(bit_data.num_elements * 64) / (double)bit_data.num_ones) / (double)SUPERBLOCK_SIZE;
    // this rounds 50% data to actaully reflect 50% data i.e 50.00001 to 50% we dont need this in interleaved because the
    // smaller superblock 63488 < 2^16 accounts for it
    log_hl_select_width = (64 - __builtin_clzll((uint64_t)((1 / float_ratio) * 0.99))) - 1;
    uint64_t total_size = build_ni_rank(bit_data);
    
    if (argv[3][0] == 's' || strcmp(argv[3], "qs") == 0) {
        total_size += build_ni_2L_select(bit_data);
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
            printf("Rank not supported: use ./ni-spider instead\n");
            exit(-1);
            break;
        case 's':
            fp_speed_test(&select_hl, WARMUP_QUERIES, bit_data.num_ones, 1);
            break;
        case 'q':

            if (argv[3][1] == 'r') {
                printf("Rank not supported: use ./ni-spider instead\n");
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
    free(hl_select_a);
    free(ll_select_a);
    free(l0_a);
    free(l1_a);


}
//!SECTION
