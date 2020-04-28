#ifndef READNIMSBIN_H
#define READNIMSBIN_H 1

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <libmseed.h>

int get_chan_name (float freq, int chan_index, char *chan_name);
int read_bin_file (FILE *f, float *s_rate, int *nscans, int *start_time,
		   int32_t **data, int *missingdataflag);

#ifdef __cplusplus
}
#endif

#endif /* READNIMSBIN_H */
