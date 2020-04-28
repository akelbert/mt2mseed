/***************************************************************************
 * readNIMSbin.c
 *
 * Functions to read the magnetotelluric channel data from the binary files
 * as output by John Booker's nimsread program. Also includes functions to
 * compute SEED standard channel names.
 *
 * Written by Anna Kelbert, COAS/OSU
 * 
 * Modidied from getinfo.c by Lana Erofeev, COAS/OSU
 *
 * last modified 8 July 2008, Chad Trabant
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libmseed.h>

/* ======================================================================= */
int get_chan_name (float freq, int chan_index, char *chan_name)
{
	//  1 HZ example:
	//	char chan_name[5][4]={{"LFN\0"},
	//	                {"LFE\0"},
	//	                {"LFZ\0"},
	//	                {"LQN\0"},
	//	                {"LQE\0"}
	//	               };
	
	char BandCode;

	if (10 <= freq && freq < 80) {
		BandCode = 'B';
	} else if (1.01 <= freq && freq < 10) {
		BandCode = 'M';
	} else if (0.5 <= freq && freq < 1.01) {
		BandCode = 'L';
	} else if (0.05 <= freq && freq < 0.5) {
		BandCode = 'V';
	} else if (0.001 <= freq && freq < 0.05) {
		BandCode = 'U';
	} else if (1e-4 <= freq && freq < 1e-3) {
		BandCode = 'R';
	} else if (1e-5 <= freq && freq < 1e-4) {
		BandCode = 'P';
	} else if (1e-5 <= freq && freq < 1e-5) {
		BandCode = 'T';
	} else if (freq < 1e-6) {
		BandCode = 'Q';
	} else {
		printf("Please refer to Appendix A of SEED manual for correct band code");
		return 0;
	}

	/* adding a "\0" to stop the string. I use pointer to access the array. 
	 the string will be the while array if the "\0" is not there.
	*/
	if (chan_index == 1) {
		sprintf(chan_name,"%cFN",BandCode);
	} else if (chan_index == 2) {
		sprintf(chan_name,"%cFE",BandCode);
	} else if (chan_index == 3) {
		sprintf(chan_name,"%cFZ",BandCode);
	} else if (chan_index == 4) {
		sprintf(chan_name,"%cQN",BandCode);
	} else if (chan_index == 5) {
		sprintf(chan_name,"%cQE",BandCode);
	} else {
		printf("We only know channel names for the first five channels");
		return 0;		
	}

	return 1;
}

/* ======================================================================= */
int read_bin_file (FILE *f, float *s_rate, int *nscans, int *start_time, 
		   int32_t **data, int *missingdataflag)
{
	/* Reads from the Fortran binary nimsread output *.bin file.
	 * Note that both Fortran and Matlab programs that write those files
	 * use 32-bit floating point and integers, irrespective of the platform
	 * on which the file was written. Therefore, instead of using
	 * sizeof(long int) for 32-bit platforms and sizeof(int) for 64-bit,
	 * use 4 for the length of an integer record in all fread statements */

	/* pointers and record length indicators */
    int32_t *itmp;
	int32_t rl, padded_rl, j1, j2;
	/* meta data from header */
	float lat, lon, decl, dt, elev;
	//int start_time[6];
	int32_t clock_zero[6];
	int32_t gaptyp, ngaps, nskip;
	int32_t *gaps;
	int swapflag = 0;
	int idx;
	/* main data block */
	//long int *data;
	
	/* read the header length record */
	if (fread(&rl, 4, 1, f) != 1) {
		printf("ERROR reading the header length record in read_bin_file\n");
		return 0;
	}
	
	/* Check if byte swapping is needed */
	padded_rl = 5108; /* (256*5-3)*4 */
	ngaps = (rl/4 - 21)/3;
	
	if (( ngaps < 0 || ngaps > 100 ) && (rl != padded_rl)) {
	  ms_gswap4a (&rl);
	  
	  ngaps = (rl/4 - 21)/3;
	  
	  if (( ngaps < 0 || ngaps > 100 ) && (rl != padded_rl)) {
	    printf("ERROR header length invalid: %d or number of gaps %d > 100\n", rl, ngaps);
	    return 0;
	  }
	  
	  printf("Byte swapping needed\n");
	  swapflag = 1;
	}
	
	j1 = rl;
	
	/* read the lat, lon, decl, dt, elev */
	if (fread(&lat, sizeof(float), 1, f) != 1) {
		printf("ERROR reading the latitude in read_bin_file\n");
		return 0;
	}
	if ( swapflag )  ms_gswap4a (&lat);
	if (fread(&lon, sizeof(float), 1, f) != 1) {
		printf("ERROR reading the longitude in read_bin_file\n");
		return 0;
	}
	if ( swapflag )  ms_gswap4a (&lon);
	if (fread(&decl, sizeof(float), 1, f) != 1) {
		printf("ERROR reading the declination in read_bin_file\n");
		return 0;
	}
	if ( swapflag )  ms_gswap4a (&decl);
	if (fread(&dt, sizeof(float), 1, f) != 1) {
		printf("ERROR reading the sampling time in read_bin_file\n");
		return 0;
	}
	if ( swapflag )  ms_gswap4a (&dt);
	if (fread(&elev, sizeof(float), 1, f) != 1) {
		printf("ERROR reading the elevation in read_bin_file\n");
		return 0;
	}
	if ( swapflag )  ms_gswap4a (&elev);
	printf("Sampling rate: %.3f Hz\n",1/dt);
	printf("Site location: (%.3f, %.3f, %.3f)\n",lat,lon,elev);

	/* read the start time and clock zero time */
	if (fread(start_time, 4, 6, f) != 6) {
		printf("ERROR reading the start time in read_bin_file\n");
		return 0;
	}
	if ( swapflag ) { for ( idx=0; idx < 6; idx++) { ms_gswap4a(start_time+idx); } }
	if (fread(clock_zero, 4, 6, f) != 6) {
		printf("ERROR reading the clock zero time in read_bin_file\n");
		return 0;
	}
	if ( swapflag ) { for ( idx=0; idx < 6; idx++) { ms_gswap4a(clock_zero+idx); } }
	printf("Time series start time: %d-%02d-%02d %d:%d:%d\n",
			start_time[0],start_time[1],start_time[2],start_time[3],start_time[4],start_time[5]);

	/* read the number of data scans and gap information */
	if (fread(nscans, 4, 1, f) != 1) {
		printf("ERROR reading the number of data scans in read_bin_file\n");
		return 0;
	}
	if ( swapflag )  ms_gswap4a (nscans);
	if (fread(&gaptyp, 4, 1, f) != 1) {
		printf("ERROR reading the gap type in read_bin_file\n");
		return 0;
	}
	if ( swapflag )  ms_gswap4a (&gaptyp);
	if (gaptyp != 2005) {
		printf("WARNING - the gap type in the file is %d, but we are assuming the gaps are filled\n",gaptyp);		
	}
	if (fread(missingdataflag, 4, 1, f) != 1) {
		printf("ERROR reading the missing data flag in read_bin_file\n");
		return 0;
	}
	if ( swapflag )  ms_gswap4a (missingdataflag);
	if (fread(&ngaps, 4, 1, f) != 1) {
		printf("ERROR reading the number of gaps in read_bin_file\n");
		return 0;
	}
	if ( swapflag )  ms_gswap4a (&ngaps);
	printf("The number of gaps in the bin file is %d\n",ngaps);

	if ( ! (gaps = (int32_t *) malloc(3*ngaps*4)) ) {
	        printf("ERROR allocating memory\n");
		return 0;
	}
	if (fread(gaps, 4, 3*ngaps, f) != 3*ngaps) {
		printf("ERROR reading the gap information in read_bin_file\n");
		return 0;
	}
	if ( swapflag ) { for ( idx=0; idx < 3*ngaps; idx++) { ms_gswap4a(gaps+idx); } }

	/* read the padding of the header if any */
	nskip = 1 + rl/4 - 22 - 3*ngaps;
	if ( ! (itmp = (int *) malloc(nskip*4)) ) {
	        printf("ERROR allocating memory\n");
		return 0;
	}
	if (fread(itmp, 4, nskip, f) != nskip) {
		printf("ERROR reading the padding of the header in read_bin_file\n");
		return 0;
	}
	if (fread(&j2, 4, 1, f) != 1) {
		printf("ERROR reading the end of header record in read_bin_file\n");
		return 0;
	}
	if ( swapflag )  ms_gswap4a (&j2);
	if ( j1 != j2 ) {
		printf("ERROR reading the end of header record marker in read_bin_file\n");
		return 0;
	}

	/* read the data length record and the data */
	if (fread(&rl, 4, 1, f) != 1) {
		printf("ERROR reading the data length record in read_bin_file\n");
		return 0;
	}
	if ( swapflag )  ms_gswap4a (&rl);
	j1 = rl;
	if ( ! (*data = (int32_t *) malloc(rl)) ) {
	        printf("ERROR allocating memory\n");
		return 0;
	}
	if (fread(*data, 4, rl/4, f) != rl/4) {
	        printf("ERROR reading the data in read_bin_file\n");
		return 0;
	}
	if ( swapflag ) { for ( idx=0; idx < rl/4; idx++) { ms_gswap4a((*data)+idx); } }
	if (fread(&j2, 4, 1, f) != 1) {
		printf("ERROR reading the end of data record in read_bin_file\n");
		return 0;
	}
	if ( swapflag )  ms_gswap4a (&j2);
	if ( j1 != j2 ) {
		printf("ERROR reading the end of data record marker in read_bin_file\n");
		return 0;
	}

	*s_rate = 1/dt;
	
	return rl;
}
