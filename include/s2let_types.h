// S2LET package
// Copyright (C) 2012-2104
// Boris Leistedt & Jason McEwen

#ifndef S2LET_TYPES
#define S2LET_TYPES

typedef enum {
    /*! Use the standard normalisation of wavelets such that they fulfil the
     *  the identity relation */
    S2LET_WAV_NORM_DEFAULT,
    /*! Normalise the wavelets as if they were spin-lowered. */
    S2LET_WAV_NORM_SPIN_LOWERED,
    /*!
     * "guard" value that equals the number of usable enum values.
     * useful in loops, for instance.
     */
    S2LET_WAV_NORM_SIZE
} s2let_wav_norm_t;

#endif
