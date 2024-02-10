// use aligned malloc
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include "helper.h"

#include <time.h>
#include <stdlib.h>
#include <x86intrin.h>

// #define WARMUP_QUERIES 100000000


// 512 = 0
// 1024 = 1
// 2048 = 2
// 4096 = 3
// #define PARAM 3
// #define VERBOSE 0

// #define MARK_SIZE 16
// #define LOG_ONES_PER_SLOT (9 + PARAM)
// #define ONES_PER_SLOT (1lu << LOG_ONES_PER_SLOT)
// #define SUPERBLOCK_SIZE 63488
// #define MARK_MASK 0xffff000000000000lu
// #define BIT_MASK ~MARK_MASK

uint64_t* sparce_select_array;

// uint64_t *modified_bit_array;
// uint64_t *l0_a;
// uint64_t *select_a;
// double float_ratio;
// uint64_t log_hl_select_width;


//SECTION - BUILD SELECT 
int build_sparce_select_from_modified(bit_meta data) {
        assert(data.num_ones > 0);

        uint64_t length = (((data.num_ones - 1) / ones_per_slot) + 1) + 1; // plus one is for the dummy

        sparce_select_array = (uint64_t*)aligned_alloc(512, length * sizeof(uint64_t));

        uint64_t position = 0;
        uint64_t sma_index = 0;
        uint64_t ones_in_sma_block = 0;
        uint64_t first_one_found = 0;
        uint64_t next = 0;

        for (uint64_t i = 0; i < data.num_elements; i++) {
            uint64_t curr_uint = modified_bit_array[i];
            uint64_t update = 64;

            if ((i % 8) == 0) {
                // this is a modified thing
                curr_uint <<= 16;
                update = 48;
            }

            if (!first_one_found) {

                if (!curr_uint) { // all 0s
                    position += update;

                    continue;
                } else {
                    // find the position of the first one in the sma array and save its position

                    uint16_t block_select = pdep_select64(curr_uint, 1);

                    sparce_select_array[sma_index] = position + block_select;
                    sma_index++;

                    // counting the rest of ones in the block, and updating the
                    // position counter to reflect that block i has been counted
                    // and setting the first_one_found flag to true
                    ones_in_sma_block += __builtin_popcountll(curr_uint) - 1;

                    position += update;
                    first_one_found = 1;
                    continue;
                }
            }

            next = __builtin_popcountll(curr_uint);


            if (next + ones_in_sma_block < ones_per_slot) {
            // position += 64;
                ones_in_sma_block += next;
            } else {
                uint16_t ones_left = ones_per_slot - ones_in_sma_block;
                uint16_t block_select = pdep_select64(curr_uint, ones_left);

                sparce_select_array[sma_index] = block_select + position;
                sma_index++;
                ones_in_sma_block = next - ones_left;
                // position += 64;
            }   


            position += update;
        }

        sparce_select_array[sma_index] = data.nbits;



    return ((sma_index + 1) * 64);
}




//!SECTION

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

    uint64_t l0_guess = sparce_select_array[ones >> log_ones_per_slot] / SUPERBLOCK_SIZE;


    uint64_t select_i = (ones - 1) >> log_ones_per_slot;

    uint64_t s_curr = sparce_select_array[select_i]; 
    uint64_t s_next = sparce_select_array[select_i + 1];


    // uint64_t basic_guess = ((((ones - ((select_i * ones_per_slot))) * (s_next - s_curr)) >> log_ones_per_slot) + 248) / 496;
    // uint64_t curr512_i = ((s_curr) / 496) + (basic_guess);

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
    clock_t b_start = clock();
    bit_meta pack = modify_package_byte_file(argv[1], atoll(argv[2]), 6);
    uint64_t total_size = (pack.num_elements * 64) - pack.nbits + (pack.l0_size * 64);


    if (argv[3][0] == 's' || argv[3][0] == 'd') {
        // meta64 tester = safe_package_byte_file(argv[1], atoll(argv[2]));
        total_size += build_sparce_select_from_modified(pack);
        // total_size += build_agressive(tester, pack);

    
    }
    clock_t b_end = clock();
    #ifndef VERBOSE
        printf("%lf\n", (double)(b_end - b_start) / (double)CLOCKS_PER_SEC);
        printf("Percent size: %lf\n", ((double)total_size * 100) / (double)(pack.nbits));
        // printf("Run Time: ");
    #endif
	// printf("%lf\n", (double)(b_end - b_start) / (double)CLOCKS_PER_SEC);
	// printf("%lf\n", ((double)total_size * 100) / (double)(pack.nbits));

    if (argc != 4) {
        help();
    }

    switch(argv[3][0]) {
        case 'r':
            printf("not supported\n");
            break;
        case 's':
            fp_speed_test(&select_hl, WARMUP_QUERIES, pack.num_ones);
            // hlrd_select_speed_test(WARMUP_QUERIES, pack.num_ones);
            break;

        case '?':
            help();
            break;
    }

    free(modified_bit_array);
    free(sparce_select_array);
    free(l0_a);

}
//!SECTION
