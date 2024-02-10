#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <x86intrin.h>

#define WARMUP_QUERIES 100000000

#define VERBOSE 0

#define MARK_SIZE 16
#define SUPERBLOCK_SIZE 63488
#define HALF_SUPERBLOCK_SIZE 31744
#define MARK_MASK 0xffff000000000000lu
#define BIT_MASK ~MARK_MASK

typedef uint64_t (*FUNC_PTR)(uint64_t);

typedef struct {
   uint64_t num_elements;
   uint64_t num_ones;
   uint64_t nbits;
   uint64_t l0_size;
   double float_ratio;

} bit_meta; 

typedef struct {
   uint64_t* bit_array;
   uint64_t num_elements;
} correctness_meta; 

uint64_t* modified_bit_array;
uint64_t*l0_a;
uint16_t* ll_select_array;
uint64_t* hl_select_array;
uint64_t log_hl_select_width;
uint64_t ones_per_slot;
uint64_t log_ones_per_slot;

int free_metadata() {

    free(modified_bit_array);
    free(hl_select_array);
    free(ll_select_array);
    free(l0_a);
    return 0;
}

int help() {
    printf("Usage: ./spider <bitfile> <numbits> <r/s>\n");
    exit(-1);

    return 0;
}

static inline uint64_t pdep_select64(uint64_t x, unsigned n) {

    return _lzcnt_u64(_pdep_u64(1lu << (__builtin_popcountll(x) - n), x));
}

//SECTION - build rank

bit_meta modify_package_byte_file(char* filename, uint64_t file_size, int log_sma_bits) {

    uint64_t size = 0;
    FILE* byte_file;
    byte_file = fopen (filename, "rb");

    if (byte_file == NULL) {
        perror("Couldn't open file");
    }

    if ((int64_t)file_size == -1) {
        fseek(byte_file, 0L, SEEK_END);
        size = ftell(byte_file) * 8; // bytes to bits
        rewind(byte_file);

    } else {
        size = file_size;
    }

    uint64_t array_length = ((size -1) / 62 + 1);
    modified_bit_array = (uint64_t*)aligned_alloc(512, (array_length + 8) * sizeof(uint64_t));
    uint64_t num_elements = 0;
    uint64_t block;
    unsigned char byte;
    uint64_t total_count = 0;
    uint64_t relative_count = 0;
    uint64_t l0_i = 0;
    uint64_t num_bits_read = 0;
    l0_a = (uint64_t*)aligned_alloc(512, sizeof(uint64_t) * (1 + 1 + ((size - 1) / SUPERBLOCK_SIZE)));
   
    for (uint64_t j = 0; j < (((size - 1) / 62) + 1); j++) { // 62 = ((x / 64) * 512) / 496
        block = 0;
        int max = 8;
        int adjust = 0;
    
        if (j % 1024 == 0) {

            l0_a[l0_i] = total_count;
            l0_i++;
            relative_count = 0;
        } 

        // each byte makes up a 8 bit section of the block
        if (j % 8 == 0) { // 8 = 512 / 64
            block |= relative_count;
            max = 6;
            adjust = __builtin_popcountll(block);
        }

        for (int i = 0; i < max; i++) {

            block = block << 8;
            byte = fgetc(byte_file);
            num_bits_read += 8;
            if (num_bits_read <= size) {
                block |= byte;
            }
      
        }

        total_count += __builtin_popcountll(block) - adjust;
        relative_count += __builtin_popcountll(block) - adjust;
      

        modified_bit_array[num_elements] = block;
        num_elements++;
    }

    // needs to be 8 because we look 8 in the future
    for (int i = num_elements; i < num_elements + 8; i++) {
        modified_bit_array[i] = -1;
    }
    num_elements += 8;
    bit_meta metadata_to_send = {};
    uint64_t num_bits = size;
    metadata_to_send.float_ratio = ((double)(num_bits) / (double)total_count) / (double)SUPERBLOCK_SIZE;


    log_hl_select_width = (64 - __builtin_clzll((uint64_t)(1 / metadata_to_send.float_ratio)));

    int base_log_ones_per_slot = 8 + log_sma_bits;
    uint64_t base_ones_per_slot = 1 << base_log_ones_per_slot;
    double sparcity = ((double)total_count / (double)(num_bits)) * 0.99;

    uint64_t rough_ones_per_slot = (uint64_t)(base_ones_per_slot * sparcity);
    log_ones_per_slot = (64 - __builtin_clzll(rough_ones_per_slot));
    ones_per_slot = 1lu << log_ones_per_slot;
    l0_a[l0_i] = l0_a[l0_i - 1] + (uint64_t)(SUPERBLOCK_SIZE / (metadata_to_send.float_ratio * SUPERBLOCK_SIZE));

    metadata_to_send.l0_size = l0_i;
    metadata_to_send.num_ones = total_count;
    // this is number of 64 bit uints blocks
    metadata_to_send.num_elements = num_elements;
    metadata_to_send.nbits = size;

    return metadata_to_send;
}

//!SECTION

//SECTION - BUILD SELECT 
int build_select_from_modified(bit_meta data) {
        assert(data.num_ones > 0);
        uint64_t hl_length = (((data.num_ones - 1) >> log_hl_select_width) + 1) + 1; // plus one is for the dummy
        uint64_t ll_length = (((data.num_ones - 1) / ones_per_slot) + 1) + 1;
        hl_select_array = (uint64_t*)aligned_alloc(512, hl_length * sizeof(uint64_t));
        ll_select_array = (uint16_t *)aligned_alloc(512, ll_length * sizeof(uint16_t));

        uint64_t hl_position = 0;
        uint64_t ll_position = 0;

        uint64_t hl_sma_index = 0;
        uint64_t ll_sma_index = 0;
        uint64_t ll_ones_in_sma_block = 0;
        uint64_t hl_ones_in_sma_block = 0;
        uint64_t first_one_found = 0;
        uint64_t next = 0;

        for (uint64_t i = 0; i < data.num_elements; i++) {
            uint64_t curr_uint = modified_bit_array[i];
            uint64_t update = 64;

            if ((i % 8) == 0) {
                // this is a modified block
                curr_uint <<= 16;
                update = 48;
            }
            
            if (i % ((1024)) == 0) {
            // allocate space for the next sma and zero out the position indicator
            ll_position = 0;
        }
            if (!first_one_found) {
            
                if (!curr_uint) { // all 0s
                    hl_position += update;
                    ll_position += update;

                    continue;
                } else {

                    // find the position of the first one in the sma array and save its position
                    uint16_t block_select = pdep_select64(curr_uint, 1);

                    hl_select_array[hl_sma_index] = (hl_position + block_select + HALF_SUPERBLOCK_SIZE) / SUPERBLOCK_SIZE;
                    ll_select_array[ll_sma_index] = ll_position + block_select;
                    hl_sma_index++;
                    ll_sma_index++;

                    // counting the rest of ones in the block, and updating the
                    // position counter to reflect that block i has been counted
                    // and setting the first_one_found flag to true
                    ll_ones_in_sma_block += __builtin_popcountll(curr_uint) - 1;
                    hl_ones_in_sma_block += __builtin_popcountll(curr_uint) - 1;

                    hl_position += update;
                    ll_position += update;
                    first_one_found = 1;
                    continue;
                }
            }

            next = __builtin_popcountll(curr_uint);
            if (next + hl_ones_in_sma_block < (1lu << log_hl_select_width)) {
                hl_ones_in_sma_block += next;
            } else {
                uint64_t hl_ones_left = (1lu << log_hl_select_width) - hl_ones_in_sma_block;
                uint64_t hl_block_select = pdep_select64(curr_uint, hl_ones_left);

                hl_select_array[hl_sma_index] = (hl_block_select + hl_position + HALF_SUPERBLOCK_SIZE) / SUPERBLOCK_SIZE;
                hl_sma_index++;
                hl_ones_in_sma_block = next - hl_ones_left;
            }

            if (next + ll_ones_in_sma_block < ones_per_slot) {
                ll_ones_in_sma_block += next;
            } else {
                uint16_t ones_left = ones_per_slot - ll_ones_in_sma_block;
                uint16_t block_select = pdep_select64(curr_uint, ones_left);

                ll_select_array[ll_sma_index] = block_select + ll_position;
                ll_sma_index++;
                ll_ones_in_sma_block = next - ones_left;
            }   


            hl_position += update;
            ll_position += update;
        }

    hl_select_array[hl_sma_index] = data.l0_size - 1;

    uint16_t ll_dummy = ll_select_array[ll_sma_index - 1];
    if ((uint64_t)(ll_dummy) + ((uint64_t)(data.float_ratio * SUPERBLOCK_SIZE) * ones_per_slot) < SUPERBLOCK_SIZE) {
        ll_select_array[ll_sma_index] = ll_dummy + ((uint64_t)(data.float_ratio * SUPERBLOCK_SIZE) * ones_per_slot);
    
    } else {
        ll_select_array[ll_sma_index] = 0;
    }

    return (hl_sma_index * 64) + (ll_sma_index * 16);
    
}

//!SECTION

//SECTION - Correctness
uint64_t count_ones(correctness_meta bit_data) {

    uint64_t count = 0;
    for (int i = 0; i < bit_data.num_elements; i++) {
        count += (uint64_t) __builtin_popcountll(bit_data.bit_array[i]);
    }
    return count;
}

// -1 find the file size manually
// this build is only for running correctness tests and does not relate to the actual structure
correctness_meta package_byte_file_for_correctness(char* filename, uint64_t file_size) {
   FILE* byte_file;
   uint64_t size = 0;

   byte_file = fopen (filename, "rb");

   if (byte_file == NULL) {
      perror("Couldn't open file");
   }


   if ((int64_t)file_size == -1) {
      fseek(byte_file, 0L, SEEK_END);
      size = ftell(byte_file);
      size = size * 8;
      rewind(byte_file);
   } else {
      size = file_size;
   }


   correctness_meta metadata_to_send = {};

   //allocate memory for the metadata array
    uint64_t* metadata_array = (uint64_t*)malloc((size / 64) * sizeof(uint64_t));

   uint64_t num_elements = 0;
   uint64_t block;
   unsigned char byte;


   for (uint64_t j = 0; j < ((size - 1) / 64) + 1; j++) {
      block = 0;

      // each byte makes up a 8 bit section of the block
      for (int i = 0; i < 8; i++) {

         block = block << 8;
         byte = fgetc(byte_file);
         block |= byte;
         
      }
      
      metadata_array[num_elements] = block;
      num_elements++;
   }

  
   metadata_to_send.bit_array = metadata_array;
   metadata_to_send.num_elements = num_elements;

   uint64_t total_count = count_ones(metadata_to_send);


   return metadata_to_send;
}

uint64_t fp_correct_test(FUNC_PTR fun, correctness_meta byte_file, uint64_t end, int is_rank_test) {
   uint64_t count = 0;
   uint64_t position = 0;
   
   for (uint64_t i = 0; i < byte_file.num_elements; i++) {
      for (int j = 63; j >= 0; j--) {
        if (is_rank_test) {
            count += !!(byte_file.bit_array[i] & (1lu << j));

            if (count != fun(position)) {
                printf("Query: %lu\n", position);
                printf("Test: %lu\n", count);
                printf("User Rank: %lu\n", fun(position));
                printf("FAIL\n");
                exit(-1);
            }

        } else {
            if (!!(byte_file.bit_array[i] & (1lu << j))) {
                    count++;
                


                if (position != fun(count)) {
                    printf("Query: %lu\n", count);
                    printf("Correct: %lu\n", position);
                    printf("User Select: %lu\n", fun(count));
                    printf("FAIL\n");
                    exit(-1);
                }

            }
        }
        position++;
        end--;
        if (!end) {
            printf("All Queries are Correct\n");

            return 0;
        }
      }
   }
    printf("All Queries are Correct\n");

   return 1;
}



//!SECTION


//SECTION - Speed Tests

// mod is num_ones for select and num_bits for rank
int fp_speed_test(FUNC_PTR fun, uint64_t num_queries, uint64_t mod) {
    
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
        numbers[i] = my_rand + 1;
    }

    for (int i = 0; i < num_queries; i++) {
        uint64_t my_rand = rand();
        my_rand = ((my_rand << 32) | rand()) % (mod);
        numbers1[i] = my_rand + 1;
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
    avg_time = (((((double)clock_time) * 1000000000)) / num_queries) / CLOCKS_PER_SEC;
    printf("%lf\n", avg_time);
    free(numbers1);

    return 0;
}

//!SECTION