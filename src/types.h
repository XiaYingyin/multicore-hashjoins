/**
 * @file    types.h
 * @author  Cagri Balkesen <cagri.balkesen@inf.ethz.ch>
 * @date    Tue May 22 16:43:30 2012
 * @version $Id: types.h 3017 2012-12-07 10:56:20Z bcagri $
 * 
 * @brief  Provides general type definitions used by all join algorithms.
 * 
 * 
 */
#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>

/**
 * @defgroup Types Common Types
 * Common type definitions used by all join implementations.
 * @{
 */

#define DIM_NUM 4

#ifdef KEY_8B /* 64-bit key/value, 16B tuples */
typedef int64_t intkey_t;
typedef int64_t value_t;
#else /* 32-bit key/value, 8B tuples */
typedef int32_t intkey_t;
typedef int32_t intvector_t;   //--define predicate vector type and change vectore type by zys
#endif
typedef int32_t value_t;
typedef int8_t vectorkey_t;
typedef struct tuple_t    tuple_t;
typedef struct column_t   column_t;
typedef struct vector_t   vector_t;
typedef struct vector_para vector_para;
typedef struct relation_t relation_t;

/** Type definition for a tuple, depending on KEY_8B a tuple can be 16B or 8B */
struct tuple_t {
    intkey_t key;
    value_t  payload;
};


/** Type definition for a column store tuple, value represents the key inside the column by zys */
struct column_t {
    intkey_t * column;
    value_t *payload;
    int8_t *bitmap;
    uint32_t  num_tuples;
};

struct vector_t {
    vectorkey_t * column;
};
struct vector_para{
    double selectivity;
    int32_t num_groups;
    uint32_t  num_tuples;
    int8_t   bitmapflag;
};

/**
 * Type definition for a relation. 
 * It consists of an array of tuples and a size of the relation.
 */
struct relation_t {
  tuple_t * tuples;
  uint32_t  num_tuples;
};

typedef struct relation_nsm_t {
	tuple_t **tuples;
	uint32_t num_columns;
	uint32_t *num_tuples;
	uint8_t **bitmaps;
	tuple_t *mem;
} relation_nsm_t;

typedef struct relation_dsm_t {
	column_t *columns;
	uint32_t num_columns;
} relation_dsm_t;

typedef int64_t (*JoinFunction)(const relation_t * const,
                                const relation_t * const,
                                relation_t * const);

typedef struct thd_param {
	relation_t * relR;
	relation_t * relS;
	JoinFunction jf;
	int nthreads;
	int result;
} thd_param_t;
/** @} */

#endif /* TYPES_H */
