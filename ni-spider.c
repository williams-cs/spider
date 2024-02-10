// use aligned malloc
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

    l0_a[l0_index] = data.num_ones;
    l1_a[l1_index] = -1;



    return ((l0_index + 1) * 64) + ((l1_index + 1) * 16);
}

//SECTION - BUILD SELECT 
int build_select(meta64 data) {
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
    // WARMUP
    for (int i = 0; i < WARMUP_QUERIES; i++) 
        result = rank(numbers[i]);
    

    start = clock();
    // TEST
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
    uint64_t count = 0;
    uint64_t position = 0;

    // iterate through every uint in file:
    for (int i = 0; i < byte_file.num_elements; i++) {
        // count every 1 seen in current uint
        for (int j = 63; j >= 0; j--) {
            count += !!(byte_file.bit_array[i] & (1lu << j));

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
                printf("All Queries are Correct\n");
                return 0;
            }
        }
    }
    printf("All Queries are Correct\n");
    return 0;
}


uint64_t cumm_select_test_hl(meta64 byte_file, uint64_t end) {
   uint64_t count = 0;
   uint64_t position = 0;
   
   for (uint64_t i = 0; i < byte_file.num_elements; i++) {
      for (int j = 63; j >= 0; j--) {
        if (!!(byte_file.bit_array[i] & (1lu << j))) {
            count++;

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
            printf("All Queries are Correct\n");

            return 0;
        }
      }
   }
   printf("All Queries are Correct\n");
   return 1;
}

//!SECTION


int help() {
    printf("Usage: ./ni-spider <bitfile> <numbits> <r/s>\n");
    exit(-1);
}

//SECTION - Main
int main(int argc, char *argv[]) {
    uint64_t query;
    clock_t b_start = clock();
    meta64 bit_data = safe_package_byte_file(argv[1], atoll(argv[2]));
    bit_a = bit_data.bit_array;
    uint64_t total_size = build_rank(bit_data);
    
    if (argv[3][0] == 's' || argv[3][0] == 'd') {
        total_size += build_select(bit_data);
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
            rank_speed_test(WARMUP_QUERIES, bit_data.nbits);
            break;
        case 's':
            hl_select_speed_test(WARMUP_QUERIES, bit_data.num_ones);
            break;
        case 'c':
            printf("%s\n", argv[1]);
           cumm_rank_test(bit_data, 100000000000lu);
            break;
        case 'd':
            printf("%s\n", argv[1]);
            cumm_select_test_hl(bit_data, 100000000000lu);
            
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
                build_select(bit_data);
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
        case '?':
            help();
            break;
    }

    free(bit_a);
    free(select_a);
    free(l0_a);
    free(l1_a);


}
//!SECTION
