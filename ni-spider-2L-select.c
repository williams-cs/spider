// use aligned malloc
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <x86intrin.h>

#define WARMUP_QUERIES 100000000


// #define VERBOSE 0


#define SUPERBLOCK_SIZE 65536
#define BASICBLOCK_SIZE 512


// uint64_t up = 0;
// uint64_t down = 0;

uint64_t* bit_a;
uint64_t *l0_a;
uint16_t *l1_a;

uint16_t *ll_select_a;
uint64_t* hl_select_a;
// double float_ratio;
uint64_t log_hl_select_width;
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

   //allocate memory for the metadata array
   uint64_t current_length = sizeof(uint64_t);
   uint64_t* metadata_array = (uint64_t*)malloc(current_length);

   uint64_t num_elements = 0;
   uint64_t block;
   unsigned char byte;
   // FILE* byte_file;

   // byte_file = fopen (filename, "rb");

   int test = 1; 

   for (uint64_t j = 0; j < ((size - 1) / 64) + 1; j++) {
      block = 0;

      // each byte makes up a 8 bit section of the block
      for (int i = 0; i < 8; i++) {

         block = block << 8;
         byte = fgetc(byte_file);
         if(test) {
            // printf("%x\n", byte);
            test = 0;
         }
         block |= byte;
         
      }
      // cast to int64_t so it can be signed (doesn't change the bits)
      // if ((int64_t)block == EOF) {
      //    break;
      // }
      
      metadata_array[num_elements] = block;
      num_elements++;

      // doubles the length of the metadata array, allocates more space to it
      if(num_elements == current_length / sizeof(uint64_t)) {
         current_length *= 2;
         
         // we want buffer in case of realloc errors
         uint64_t* buffer = (uint64_t*)realloc(metadata_array, current_length);
         if (buffer == NULL) {

            perror("realloc error");
         }
         metadata_array = buffer;
      }
   }

   metadata_to_send.nbits = size;
   metadata_to_send.bit_array = metadata_array;
   metadata_to_send.num_elements = num_elements;

   uint64_t total_count = count_ones(metadata_to_send);
   
    int base_log_ones_per_slot = 12;
    uint64_t base_ones_per_slot = 1 << base_log_ones_per_slot;
    uint64_t rough_ones_per_slot = (uint64_t)(base_ones_per_slot * ((double)total_count / (double)(num_elements * 64)) * .99);
    // printf("ratL %ld\n", (uint64_t)((double)(num_elements * 64) / (double)total_count ));
    log_ones_per_slot = (64 - __builtin_clzll(rough_ones_per_slot));
    double sparcity = (double)total_count / (double)(num_elements * 64) * 0.99;
    // printf("Sparcity: %lf\n", sparcity);
    // printf("log ones per slot: %ld\n", log_ones_per_slot);
    ones_per_slot = 1lu << log_ones_per_slot;


   metadata_to_send.float_ratio = ((double)(num_elements * 64) / (double)total_count) / (double)SUPERBLOCK_SIZE;
   // printf("%lf\n", metadata_to_send.float_ratio);
   

   metadata_to_send.num_ones =total_count;
   


   return metadata_to_send;
}


static inline uint64_t pdep_select64(uint64_t x, unsigned n) {
    // printf("%lld\n",_lzcnt_u64(_pdep_u64(1lu << (__builtin_popcountll(x) - n), x)));

    return _lzcnt_u64(_pdep_u64(1lu << (__builtin_popcountll(x) - n), x));
}

uint64_t build_rank(meta64 data) {
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

    l0_a[l0_index] = data.nbits;
    l1_a[l1_index] = -1;



    return ((l0_index + 1) * 64) + ((l1_index + 1) * 16);
}

//SECTION - BUILD SELECT 
uint64_t build_select(meta64 data) {
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
            // ones_in_sma_block = 0;
            // first_one_found = 0;
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
                // position += 64;
                hl_ones_in_sma_block += next;
            } else {
                uint64_t hl_ones_left = (1lu << log_hl_select_width) - hl_ones_in_sma_block;
                uint64_t hl_block_select = pdep_select64(curr_uint, hl_ones_left);

                hl_select_a[hl_sma_index] = hl_block_select + hl_position;
                hl_sma_index++;
                hl_ones_in_sma_block = next - hl_ones_left;
            // position += 64;
            }

            if (next + ll_ones_in_sma_block < ones_per_slot) {
            // position += 64;
                ll_ones_in_sma_block += next;
            } else {
                uint16_t ones_left = ones_per_slot - ll_ones_in_sma_block;
                uint16_t block_select = pdep_select64(curr_uint, ones_left);

                ll_select_a[ll_sma_index] = block_select + ll_position;
                ll_sma_index++;
                ll_ones_in_sma_block = next - ones_left;
                // position += 64;
            }   


            hl_position += update;
            ll_position += update;
        }

    uint64_t hl_dummy = hl_select_a[hl_sma_index - 1];// + (total_ratio * ONES_PER_SLOT) / 512;
    hl_select_a[hl_sma_index] = hl_dummy + ((uint64_t)(data.float_ratio * SUPERBLOCK_SIZE) << log_hl_select_width);

    uint16_t ll_dummy = ll_select_a[ll_sma_index - 1];// + (total_ratio * ONES_PER_SLOT) / 512;
    if ((uint64_t)(ll_dummy) + ((uint64_t)(data.float_ratio * SUPERBLOCK_SIZE) * ones_per_slot) < (1 << 16)) {
        ll_select_a[ll_sma_index] = ll_dummy + ((uint64_t)(data.float_ratio * SUPERBLOCK_SIZE) * ones_per_slot);
    
    } else {
        ll_select_a[ll_sma_index] = 0;
    }

    return ((hl_sma_index + 1) * 64) + ((ll_sma_index + 1) * 16);
}




//!SECTION
// this looks like poppy copy paste (change)
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
// uint64_t rank(uint64_t pos) {
//     uint64_t section = pos / 496;
//     uint64_t l0_i = section >> 7;
//     uint64_t start = section << 3;

//     return l0_a[l0_i] + ((mod_a[start] & MARK_MASK) >> 48) + rank_512(mod_a, start, (pos - (section * 496)) + 16);
// }

uint64_t rank(uint64_t pos) {
    uint64_t l0_i = pos >> 16;
    uint64_t l1_i = pos >> 9;
    uint64_t start = l1_i << 3;
    uint64_t so_far = l0_a[l0_i] + l1_a[l1_i];
    // printf("%lu\n", pos - so_far);

    return so_far + rank_512(bit_a, start, pos - (l1_i * BASICBLOCK_SIZE));
}


uint64_t select_hl(uint64_t ones) {

    uint64_t l0_guess = (hl_select_a[(ones >> log_hl_select_width)] + 32768) / SUPERBLOCK_SIZE;
    // printf("l0g: %lu\n", l0_guess);

    while (ones <= l0_a[l0_guess]) {
        // up++;
        l0_guess--;
    }

    while (ones > l0_a[l0_guess + 1]) {
        // down++;

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

    // uint64_t basic_guess = ((((ones - ((select_i * ones_per_slot))) * (s_next - s_curr)) >> log_ones_per_slot) + 256) / 512;
    // uint64_t curr512_i = (((s_curr + l0_guess * SUPERBLOCK_SIZE)) / 512) + (basic_guess);


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


//SECTION - Speed Tests
int rank_speed_test(uint64_t num_queries, uint64_t nbits) {
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
        my_rand = ((my_rand << 32) | rand()) % (nbits);
        numbers[i] = my_rand;
    }


    for (int i = 0; i < num_queries; i++) {
        uint64_t my_rand = rand();
        my_rand = ((my_rand << 32) | rand()) % (nbits);
        numbers1[i] = my_rand;
    }
    volatile uint64_t result;
    // warmup'
    for (int i = 0; i < WARMUP_QUERIES; i++) 
        result = rank(numbers[i]);
    

    start = clock();

    for (int i = 0; i < num_queries; i++) 
        result = rank(numbers1[i]);

    end = clock();

    clock_t clock_time = end - start;
    avg_time = (((((double)clock_time) * 1000000000)) / num_queries) / CLOCKS_PER_SEC;
    printf("%lf\n", avg_time);
    free(numbers1);
    free(numbers);
    return 1;
}

int hl_select_speed_test(uint64_t num_queries, uint64_t num_ones) {
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
        my_rand = ((my_rand << 32) | rand()) % (num_ones);
        numbers[i] = my_rand + 1;
    }

    for (int i = 0; i < num_queries; i++) {
        uint64_t my_rand = rand();
        my_rand = ((my_rand << 32) | rand()) % (num_ones);
        numbers1[i] = my_rand + 1;
    }

    volatile uint64_t result;
    // WARMUP
    for (int i = 0; i < WARMUP_QUERIES; i++)
        result = select_hl(numbers[i]);


    // TEST
    start = clock();
    for (int i = 0; i < num_queries; i++)
        result = select_hl(numbers1[i]);
    

    end = clock();
    clock_t clock_time = end - start;
    avg_time = (((((double)clock_time) * 1000000000)) / num_queries) / CLOCKS_PER_SEC;
    printf("%lf\n", avg_time);
    free(numbers1);
    free(numbers);
    return 1;
}

//!SECTION

//SECTION - Correctness tests
// cumulative test for rank: performs rank queries from 0 to end
int cumm_rank_test(meta64 byte_file, uint64_t end) {
    // count of every 1 seen up to current index (stored in position)
    uint64_t count = 0;
    // current index
    uint64_t position = 0;

    // iterate through every uint in file:
    for (int i = 0; i < byte_file.num_elements; i++) {
        // count every 1 seen in current uint
        for (int j = 63; j >= 0; j--) {
            count += !!(byte_file.bit_array[i] & (1lu << j));

            // if (position % 10000 == 0) {
            //     printf("Query: %lu\n", position);
            // }
            // check that manual count is equal to our rank query program
            if (count != rank(position)) {
                printf("Query: %lu\n", position);
                printf("Test: %lu\n", count);
                printf("ours: %lu\n", rank(position));
                printf("FAIL\n");
                exit(-1);
            }
            position++;
            end--;
            if (end == 0) {
                printf("CLEAN\n");
                return 0;
            }
        }
    }
    printf("CLEAN\n");
    return 0;
}


uint64_t cumm_select_test_hl(meta64 byte_file, uint64_t end) {
   uint64_t count = 0;
   uint64_t position = 0;
   
   for (uint64_t i = 0; i < byte_file.num_elements; i++) {
      for (int j = 63; j >= 0; j--) {
        if (!!(byte_file.bit_array[i] & (1lu << j))) {
            count++;
            if (count % 10000000 == 0) {
                // printf("Query: %lu\n", count);
            }

            if (position != select_hl(count)) {
                printf("Query: %lu\n", count);
                printf("Correct: %lu\n", position);
                printf("Selectx: %lu\n", select_hl(count));
                printf("FAIL\n");
                exit(-1);
            }

        }
        position++;
        end--;
        if (!end) {
            printf("CLEAN\n");

            return 0;
        }
      }
   }
   printf("CLEAN\n");
   return 1;
}

//!SECTION


int help() {
    printf("TODO\n");
    exit(-1);
}

//SECTION - Main
int main(int argc, char *argv[]) {
    clock_t b_start = clock();
    meta64 bit_data = safe_package_byte_file(argv[1], atoll(argv[2]));
    bit_a = bit_data.bit_array;
    double float_ratio = ((double)(bit_data.num_elements * 64) / (double)bit_data.num_ones) / (double)SUPERBLOCK_SIZE;
    // this rounds 50% data to actaully reflect 50% data i.e 50.00001 to 50% we dont need this in interleaved because the
    // smaller superblock 63488 < 2^16 accounts for it
    log_hl_select_width = (64 - __builtin_clzll((uint64_t)((1 / float_ratio) * 0.99))) - 1;
    uint64_t total_size = build_rank(bit_data) + build_select(bit_data);




    clock_t b_end = clock();

    if (argv[3][0] == 'c' || argv[3][0] == 'd') {
        // tester = safe_package_byte_file(argv[1], atoll(argv[2]));
    }
    #ifndef VERBOSE
        printf("%lf\n", (double)(b_end - b_start) / (double)CLOCKS_PER_SEC);
        printf("Percent size: %lf\n", ((double)total_size * 100) / (double)(bit_data.nbits));
        // printf("Run Time: ");
    #endif
	// printf("%lf\n", (double)(b_end - b_start) / (double)CLOCKS_PER_SEC);
	// printf("%lf\n", ((double)total_size * 100) / (double)(bit_data.nbits));
    // rank_speed_test(WARMUP_QUERIES, bit_data.nbits);
            // cumm_select_test_hl(bit_data, 100000000000lu);

    // cumm_rank_test(bit_data, 100000000000lu);
    // for (int i = 0; i < 10; i++) {
    //     printf("%lu\n", hl_select_a[i]);
        
    // }
    if (argc != 4) {
        help();
    }

    switch(argv[3][0]) {
        case 'r':
            rank_speed_test(WARMUP_QUERIES, bit_data.nbits);
            break;
        case 's':
            hl_select_speed_test(WARMUP_QUERIES, bit_data.num_ones);
            // hlrd_select_speed_test(WARMUP_QUERIES, pack.num_ones);
            break;
        case 'c':
            printf("%s\n", argv[1]);
           cumm_rank_test(bit_data, 100000000000lu);
            break;
        case 'd':
            printf("%s\n", argv[1]);
            cumm_select_test_hl(bit_data, 100000000000lu);
            
            break;

        case '?':
            help();
            break;
    }
    // printf("up: %lu, down: %lu\n", up , down);

    free(bit_a);
    free(hl_select_a);
    free(ll_select_a);
    free(l0_a);
    free(l1_a);


}
//!SECTION
