#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <x86intrin.h>

#define WARMUP_QUERIES 100000000
#define SUPERBLOCK_SIZE 65536
#define BASICBLOCK_SIZE 512

#define VERBOSE

uint64_t* bit_a;
uint64_t *l0_a;
uint16_t *l1_a;

uint64_t* select_a;


uint64_t log_ones_per_slot;
uint64_t ones_per_slot;

typedef struct metadata {
   uint64_t* bit_array;
   uint64_t num_elements;
   uint64_t num_ones;
   uint64_t nbits;
   double float_ratio;
} meta64; 


uint64_t count_ones(struct metadata bit_data) {

    uint64_t count = 0;
    for (int i = 0; i < bit_data.num_elements; i++) {
        count += (uint64_t) __builtin_popcountll(bit_data.bit_array[i]);
    }
    return count;
}


// -1 find the file size manually
meta64 safe_package_byte_file(char* filename, uint64_t file_size) {
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


   meta64 metadata_to_send = {};
   
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

   metadata_to_send.nbits = size;
   metadata_to_send.bit_array = metadata_array;
   metadata_to_send.num_elements = num_elements;

   uint64_t total_count = count_ones(metadata_to_send);
   
    int base_log_ones_per_slot = 14;
    uint64_t base_ones_per_slot = 1 << base_log_ones_per_slot;
    uint64_t rough_ones_per_slot = (uint64_t)(base_ones_per_slot * ((double)total_count / (double)(num_elements * 64)) * 0.99);
    log_ones_per_slot = (64 - __builtin_clzll(rough_ones_per_slot));
    double sparcity = (double)total_count / (double)(num_elements * 64) * 0.99;
    ones_per_slot = 1lu << log_ones_per_slot;

    metadata_to_send.float_ratio = ((double)(num_elements * 64) / (double)total_count) / (double)SUPERBLOCK_SIZE;
   

   metadata_to_send.num_ones =total_count;
   


   return metadata_to_send;
}

static inline uint64_t pdep_select64(uint64_t x, unsigned n) {
    return _lzcnt_u64(_pdep_u64(1lu << (__builtin_popcountll(x) - n), x));
}

uint64_t build_ni_rank(meta64 data) {
    l0_a = (uint64_t*)aligned_alloc(512, sizeof(uint64_t) * (1 + 1 + ((data.nbits - 1) / SUPERBLOCK_SIZE))); // maybe dont need plus 1
    l1_a = (uint16_t*)aligned_alloc(512, sizeof(uint16_t) * (1 + 1 + ((data.nbits - 1) / BASICBLOCK_SIZE)));
    uint64_t l1_count = 0;
    uint64_t l1_index = 0;
    uint64_t l0_index = 0;
    uint64_t l0_count = 0;


    for (int i = 0; i < data.num_elements; i++) {

        if (i % (SUPERBLOCK_SIZE / 64) == 0) {
            l0_a[l0_index] = l0_count;
            l0_index++;
            l1_count = 0;


        }

        if (i % (BASICBLOCK_SIZE / 64) == 0) {

            assert(l1_count < (1 << 16));
            l1_a[l1_index] = l1_count;
            l1_index++;
            
        }
    
        l1_count += __builtin_popcountll(data.bit_array[i]);
        l0_count += __builtin_popcountll(data.bit_array[i]);


    }

    l0_a[l0_index] = data.num_ones;
    l1_a[l1_index] = -1;



    return ((l0_index + 1) * 64) + ((l1_index + 1) * 16);
}

//SECTION - BUILD SELECT 
int build_ni_select(meta64 data) {
        assert(data.num_ones > 0);
        uint64_t length = (((data.num_ones - 1) / ones_per_slot) + 1) + 1;
        select_a = (uint64_t*)aligned_alloc(512, length * sizeof(uint64_t));

        uint64_t position = 0;

        uint64_t sma_index = 0;
        uint64_t ones_in_sma_block = 0;
        uint64_t first_one_found = 0;
        uint64_t next = 0;

        for (uint64_t i = 0; i < data.num_elements; i++) {
            uint64_t curr_uint = data.bit_array[i];
            
            if (!first_one_found) {
            
                if (!curr_uint) { // all 0s
                    position += 64;

                    continue;
                } else {
                    // find the position of the first one in the sma array and save its position

                    uint16_t block_select = pdep_select64(curr_uint, 1);

                    select_a[sma_index] = position + block_select;
                    sma_index++;

                    // counting the rest of ones in the block, and updating the
                    // position counter to reflect that block i has been counted
                    // and setting the first_one_found flag to true
                    ones_in_sma_block += __builtin_popcountll(curr_uint) - 1;

                    position += 64;
                    first_one_found = 1;
                    continue;
                }
            }

            next = __builtin_popcountll(curr_uint);

            if (next + ones_in_sma_block < ones_per_slot) {
                ones_in_sma_block += next;
            } else {
                uint64_t ones_left = ones_per_slot - ones_in_sma_block;
                uint64_t block_select = pdep_select64(curr_uint, ones_left);

                select_a[sma_index] = block_select + position;
                sma_index++;
                ones_in_sma_block = next - ones_left;
            }   


            position += 64;
        }

    select_a[sma_index] = data.nbits;
    return ((sma_index + 1) * 64);
}

//SECTION - BUILD 2L SELECT 

uint16_t *ll_select_a;
uint64_t* hl_select_a;
uint64_t log_hl_select_width;


uint64_t build_ni_2L_select(meta64 data) {
        assert(data.num_ones > 0);
        uint64_t hl_length = (((data.num_ones - 1) >> log_hl_select_width) + 1) + 1; // plus one is for the dummy
        uint64_t ll_length = (((data.num_ones - 1) / ones_per_slot) + 1) + 1;
        hl_select_a = (uint64_t*)aligned_alloc(512, hl_length * sizeof(uint64_t));
        ll_select_a = (uint16_t *)aligned_alloc(512, ll_length * sizeof(uint16_t));

        uint64_t hl_position = 0;
        uint64_t ll_position = 0;

        uint64_t hl_sma_index = 0;
        uint64_t ll_sma_index = 0;
        uint64_t ll_ones_in_sma_block = 0;
        uint64_t hl_ones_in_sma_block = 0;
        uint64_t first_one_found = 0;
        uint64_t next = 0;

        for (uint64_t i = 0; i < data.num_elements; i++) {
            uint64_t curr_uint = data.bit_array[i];
            uint64_t update = 64;
            
            if (i % ((1024)) == 0) {
            // allocate space for the next sma and zero out all
            // the position indicators
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

                    hl_select_a[hl_sma_index] = hl_position + block_select;
                    ll_select_a[ll_sma_index] = ll_position + block_select;
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

                hl_select_a[hl_sma_index] = hl_block_select + hl_position;
                hl_sma_index++;
                hl_ones_in_sma_block = next - hl_ones_left;
            }

            if (next + ll_ones_in_sma_block < ones_per_slot) {
                ll_ones_in_sma_block += next;
            } else {
                uint16_t ones_left = ones_per_slot - ll_ones_in_sma_block;
                uint16_t block_select = pdep_select64(curr_uint, ones_left);

                ll_select_a[ll_sma_index] = block_select + ll_position;
                ll_sma_index++;
                ll_ones_in_sma_block = next - ones_left;
            }   


            hl_position += update;
            ll_position += update;
        }

    uint64_t hl_dummy = hl_select_a[hl_sma_index - 1];
    hl_select_a[hl_sma_index] = hl_dummy + ((uint64_t)(data.float_ratio * SUPERBLOCK_SIZE) << log_hl_select_width);

    uint16_t ll_dummy = ll_select_a[ll_sma_index - 1];
    if ((uint64_t)(ll_dummy) + ((uint64_t)(data.float_ratio * SUPERBLOCK_SIZE) * ones_per_slot) < (1 << 16)) {
        ll_select_a[ll_sma_index] = ll_dummy + ((uint64_t)(data.float_ratio * SUPERBLOCK_SIZE) * ones_per_slot);
    
    } else {
        ll_select_a[ll_sma_index] = 0;
    }

    return ((hl_sma_index + 1) * 64) + ((ll_sma_index + 1) * 16);
}

//!SECTION


//!SECTION


int help() {
    printf("Usage: ./ni-spider <bitfile> <numbits> <r/s>\n");
    exit(-1);
}