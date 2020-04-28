/***************************************************************************
 * mt2mseed.c
 *
 * Data conversion from magnetotelluric timeseries as "bin" files to Mini-SEED.
 *
 * Written by Chad Trabant, IRIS Data Management Center
 *
 * last modified 10 July 2008, Anna Kelbert, COAS/OSU
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>

#include <libmseed.h>

#include "readNIMSbin.h"

#define VERSION "1.1"
#define PACKAGE "mt2mseed"

struct listnode {
  char *key;
  char *data;
  struct listnode *next;
};

static int packmsr (MSRecord *msr);
static int binconvert (char *binfile);
static int parameter_proc (int argcount, char **argvec);
static char *getoptval (int argcount, char **argvec, int argopt);
static int readlistfile (char *listfile);
static void addnode (struct listnode **listroot, char *key, char *data);
static void record_handler (char *record, int reclen, void *handlerdata);
static void usage (void);

static int   verbose     = 0;
static int   packreclen  = -1;
static int   encoding    = 11;
static int   byteorder   = -1;
static int   chanfiles   = 0;
static char  srateblkt   = 0;
static char *network     = "EM";
static char *station     = 0;
static char *location    = 0;
static char *outfile     = 0;
static FILE *outfp       = 0;

/* A list of input files */
struct listnode *filelist = 0;

static int64_t packedsamples = 0;
static int64_t packedrecords = 0;

int
main (int argc, char **argv)
{
  struct listnode *flp;
  
  /* Process given parameters (command line and parameter file) */
  if (parameter_proc (argc, argv) < 0)
    return -1;
  
  /* Convert each input bin file */
  flp = filelist;
  while ( flp != 0 )
    {
      if ( verbose )
	fprintf (stderr, "Reading %s\n", flp->data);
      
      binconvert (flp->data);
      
      flp = flp->next;
    }
  
  fprintf (stderr, "Packed %" PRId64 " samples into %" PRId64 " records\n",
	   packedsamples, packedrecords);
  
  /* Close user specified output file */
  if ( outfp )
    fclose (outfp);
  
  return 0;
}  /* End of main() */


/***************************************************************************
 * packmsr:
 *
 * Pack all samples in the specified MSRecord.  If a single output
 * file has been specified it will be opened and all output will be
 * written to it, otherwise filenames will be created for each channel
 * segment and will include the start time of the segment.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
packmsr (MSRecord *msr)
{
  FILE *ofp = 0;
  char ofname[1024], timestr[20];
  int64_t trpackedsamples = 0;
  int64_t trpackedrecords = 0;
  
  if ( outfile )
    {
      /* Open user specified output file */
      if ( ! outfp )
	{
	  if ( ! (outfp = fopen (outfile, "w")) )
	    {
	      fprintf (stderr, "Error opening output file: %s\n",
		       strerror(errno));
	      return -1;
	    }
	}
      
      ofp = outfp;
    }
  	else if (! chanfiles) {
      /* Generate the output file name for all channels and segments
       * and open output file */
      if ( ! outfp )
  	{
      ms_hptime2isotimestr (msr->starttime, timestr, 0);

      snprintf (ofname, sizeof(ofname), "%s.%s.%s",
    		msr->network, msr->station, timestr);

  	  if ( ! (outfp = fopen (ofname, "w")) )
  	    {
  	      fprintf (stderr, "Error opening output file: %s\n",
  		       strerror(errno));
  	      return -1;
  	    }
  	}
        
       ofp = outfp; 
  	}
  	else	  
    {
      /* Generate the output file name for new channel and/or
       * segment and open output file */
      ms_hptime2isotimestr (msr->starttime, timestr, 0);

      snprintf (ofname, sizeof(ofname), "%s.%s.%s.%s",
		msr->network, msr->station, timestr, msr->channel);
      
      if ( ! (ofp = fopen (ofname, "w")) )
	{
	  fprintf (stderr, "Error opening output file: %s\n",
		   strerror(errno));
	  return -1;
	}
    }
  
  msr->encoding = encoding;
  
  /* Pack output data */
  trpackedrecords = msr_pack (msr, &record_handler, ofp,
			      &trpackedsamples, 1, verbose-2);
  
  if ( trpackedrecords < 0 )
    {
      fprintf (stderr, "Error packing data\n");
      return -1;
    }
  else
    {
      packedrecords += trpackedrecords;
      packedsamples += trpackedsamples;
    }
  
  /* Close file only if not the user specified file */
  if ( ofp && ofp != outfp )
    {
      fclose (ofp);
    }
  
  return 0;
}  /* End of packmsr() */


/***************************************************************************
 * binconvert:
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
binconvert (char *binfile)
{
  FILE *ifp = 0;
  MSRecord *msr = 0;
  struct blkt_1000_s Blkt1000;
  struct blkt_100_s Blkt100;
  
  int dataidx;
  int datacnt;
  int startidx;
  int channel;
  int yday;
  int nscans; 
  int datasize;
  int start_time[6];
  int missingdataflag;
  hptime_t starttime;
  int32_t *idata = 0;
  int32_t *cdata = 0;
  float samprate;
  char chan[6];
  
  /* Open input file */
  if ( (ifp = fopen (binfile, "rb")) == NULL )
    {
      fprintf (stderr, "Cannot open input file: %s (%s)\n",
	       binfile, strerror(errno));
      return -1;
    }
  
  /* Parse bin file */
  if ( ! (datasize = read_bin_file (ifp, &samprate, &nscans, start_time,
				    &idata, &missingdataflag)) )
    {
      fprintf (stderr, "[%s] Error reading input bin file\n", binfile);
      return -1;
    }
  
  if ( verbose > 0 )
    fprintf (stderr, "[%s] Missing data flag (value): %d\n", binfile, missingdataflag);
  
  if ( samprate <= 0.0 )
    {
      fprintf (stderr, "[%s] Error with sample rate\n", binfile);
      return -1;
    }
  
  if ( datasize != (nscans * 5 * sizeof(int32_t)) )
    {
      fprintf (stderr, "[%s] Unexpected data array size (%d bytes) for %d scans of 5 channels",
	       binfile, datasize, nscans);
      return -1;
    }
  
  if ( ms_md2doy (start_time[0], start_time[1], start_time[2], &yday) )
    {
      fprintf (stderr, "Error converting month and day-of-month to day-of-year\n");
      fprintf (stderr, "  Input year: %d, month: %d, day-of-month: %d\n",
	       start_time[0], start_time[1], start_time[2]);
      return -1;
    }
  
  if ( verbose >= 1 )
    {
      printf ("[%s] Start time: %d,%d,%d:%d:%d\n", binfile,
	      start_time[0], yday, start_time[3], start_time[4], start_time[5]);
      printf ("[%s] Sample rate is %.3f HZ for %d data scans\n",
	      binfile, samprate, nscans);
    }
  
  /* Initialize MSRecord */
  if ( ! (msr = msr_init(msr)) )
    {
      fprintf (stderr, "[%s] Error initializing MSRecord\n", binfile);
      return -1;
    }

  /* Allocate channel specific buffer */
  if ( ! (cdata = (int32_t *) malloc (sizeof(int32_t) * nscans)) )
    {
      fprintf (stderr, "[%s] Error allocating memory\n", binfile);
      return -1;
    }
  msr->datasamples = cdata;
  msr->sampletype = 'i';
  
  /* Set start time, sample rate and number of samples */
  starttime = ms_time2hptime (start_time[0], yday, start_time[3],
			      start_time[4], start_time[5], 0);
  msr->samprate = samprate;
    
  /* Set network, station and location */
  if ( network )
    ms_strncpclean (msr->network, network, 2);
  if ( station )
    ms_strncpclean (msr->station, station, 5);
  if ( location )
    ms_strncpclean (msr->location, location, 2);
  
  /* Add blockettes 1000 to MSRecord */
  memset (&Blkt1000, 0, sizeof(struct blkt_1000_s));
  msr_addblockette (msr, (char *) &Blkt1000,
		    sizeof(struct blkt_1001_s), 1000, 0);
  
  /* Add blockette 100 to template if requested */
  if ( srateblkt )
    {
      memset (&Blkt100, 0, sizeof(struct blkt_100_s));
      Blkt100.samprate = (float) msr->samprate;
      msr_addblockette (msr, (char *) &Blkt100,
			sizeof(struct blkt_100_s), 100, 0);
    }
  
  /* Loop over 5 channels */
  for ( channel=0; channel < 5 ; channel++ )
    {
      if ( get_chan_name (samprate, channel+1, chan) != 1 )
	{
	  fprintf (stderr, "[%s] Unable to determine channel codes for channel number %d",
		   binfile, channel+1);
	  break;
	}
      
      if ( verbose > 1 )
	fprintf (stderr, "[%s] Reading data for channel %d (%s)\n",
		 binfile, channel+1, chan);
      
      /* Set channel codes */
      ms_strncpclean (msr->channel, chan, 3);
      
      dataidx = 0;
      
      while ( dataidx < nscans )
	{
	  /* Extract data array for this channel */
	  for (datacnt=0, startidx=dataidx; dataidx < nscans; dataidx++)
	    {
	      if ( idata[(5 * dataidx) + channel] == missingdataflag ||
		   idata[(5 * dataidx) + channel] >= 2147483647 )
		{
		  dataidx++;
		  break;
		}
	      
	      cdata[datacnt++] = idata[(5 * dataidx) + channel];
	    }
	  
	  if ( datacnt > 0 )
	    {
	      if ( verbose >= 1 )
		{
		  fprintf (stderr, "[%s] %d samps @ %.6f Hz for N: '%s', S: '%s', L: '%s', C: '%s'\n",
			   binfile, datacnt, msr->samprate,
			   msr->network, msr->station,  msr->location, msr->channel);
		}
	      
	      /* Set start time and sample counts */
	      msr->starttime = starttime + ((startidx / msr->samprate) * HPTMODULUS);
	      msr->samplecnt = msr->numsamples = datacnt;
	      
	      /* Pack data into records */
	      if ( packmsr (msr) )
		{
		  fprintf (stderr, "[%s] Error packing Mini-SEED\n", binfile);
		  break;
		}
	    }
	}
    }
  
  fclose (ifp);
  
  if ( idata )
    free (idata);
  
  if ( cdata )
    free (cdata);
  
  msr->datasamples = 0;
  if ( msr )
    msr_free (&msr);
  
  return 0;
}  /* End of binconvert() */


/***************************************************************************
 * parameter_proc:
 * Process the command line parameters.
 *
 * Returns 0 on success, and -1 on failure.
 ***************************************************************************/
static int
parameter_proc (int argcount, char **argvec)
{
  int optind;

  /* Process all command line arguments */
  for (optind = 1; optind < argcount; optind++)
    {
      if (strcmp (argvec[optind], "-V") == 0)
	{
	  fprintf (stderr, "%s version: %s\n", PACKAGE, VERSION);
	  exit (0);
	}
      else if (strcmp (argvec[optind], "-h") == 0)
	{
	  usage();
	  exit (0);
	}
      else if (strncmp (argvec[optind], "-v", 2) == 0)
	{
	  verbose += strspn (&argvec[optind][1], "v");
	}
      else if (strcmp (argvec[optind], "-S") == 0)
	{
	  srateblkt = 1;
	}
      else if (strcmp (argvec[optind], "-C") == 0)
	{
	  chanfiles = 1;
	}
      else if (strcmp (argvec[optind], "-n") == 0)
	{
	  network = getoptval(argcount, argvec, optind++);
	}
      else if (strcmp (argvec[optind], "-s") == 0)
	{
	  station = getoptval(argcount, argvec, optind++);
	}
      else if (strcmp (argvec[optind], "-l") == 0)
	{
	  location = getoptval(argcount, argvec, optind++);
	}
      else if (strcmp (argvec[optind], "-r") == 0)
	{
	  packreclen = strtoul (getoptval(argcount, argvec, optind++), NULL, 10);
	}
      else if (strcmp (argvec[optind], "-e") == 0)
	{
	  encoding = strtoul (getoptval(argcount, argvec, optind++), NULL, 10);
	}
      else if (strcmp (argvec[optind], "-b") == 0)
	{
	  byteorder = strtoul (getoptval(argcount, argvec, optind++), NULL, 10);
	}
      else if (strcmp (argvec[optind], "-o") == 0)
	{
	  outfile = getoptval(argcount, argvec, optind++);
	}
      else if (strncmp (argvec[optind], "-", 1) == 0 &&
	       strlen (argvec[optind]) > 1 )
	{
	  fprintf(stderr, "Unknown option: %s\n", argvec[optind]);
	  exit (1);
	}
      else
	{
	  addnode (&filelist, NULL, argvec[optind]);
	}
    }

  /* Make sure an input files were specified */
  if ( filelist == 0 )
    {
      fprintf (stderr, "No input files were specified\n\n");
      fprintf (stderr, "%s version %s\n\n", PACKAGE, VERSION);
      fprintf (stderr, "Try %s -h for usage\n", PACKAGE);
      exit (1);
    }

  /* Report the program version */
  if ( verbose )
    fprintf (stderr, "%s version: %s\n", PACKAGE, VERSION);
  
  /* Sanity check encoding */
  if ( encoding != 3 && encoding != 10 && encoding != 11 )
    {
      fprintf (stderr, "Unsupported encoding type: %d\n", encoding);
      exit (1);
    }
  
  /* Check the input files for any list files, if any are found
   * remove them from the list and add the contained list */
  if ( filelist )
    {
      struct listnode *prevln, *ln;
      char *lfname;
      
      prevln = ln = filelist;
      while ( ln != 0 )
	{
	  lfname = ln->data;
	  
	  if ( *lfname == '@' )
	    {
	      /* Remove this node from the list */
	      if ( ln == filelist )
		filelist = ln->next;
	      else
		prevln->next = ln->next;
	      
	      /* Skip the '@' first character */
	      if ( *lfname == '@' )
		lfname++;

	      /* Read list file */
	      readlistfile (lfname);
	      
	      /* Free memory for this node */
	      if ( ln->key )
		free (ln->key);
	      free (ln->data);
	      free (ln);
	    }
	  else
	    {
	      prevln = ln;
	    }
	  
	  ln = ln->next;
	}
    }

  return 0;
}  /* End of parameter_proc() */


/***************************************************************************
 * getoptval:
 * Return the value to a command line option; checking that the value is 
 * itself not an option (starting with '-') and is not past the end of
 * the argument list.
 *
 * argcount: total arguments in argvec
 * argvec: argument list
 * argopt: index of option to process, value is expected to be at argopt+1
 *
 * Returns value on success and exits with error message on failure
 ***************************************************************************/
static char *
getoptval (int argcount, char **argvec, int argopt)
{
  if ( argvec == NULL || argvec[argopt] == NULL ) {
    fprintf (stderr, "getoptval(): NULL option requested\n");
    exit (1);
    return 0;
  }
  
  /* Special case of '-o -' usage */
  if ( (argopt+1) < argcount && strcmp (argvec[argopt], "-o") == 0 )
    if ( strcmp (argvec[argopt+1], "-") == 0 )
      return argvec[argopt+1];
  
  if ( (argopt+1) < argcount && *argvec[argopt+1] != '-' )
    return argvec[argopt+1];
  
  fprintf (stderr, "Option %s requires a value\n", argvec[argopt]);
  exit (1);
  return 0;
}  /* End of getoptval() */


/***************************************************************************
 * readlistfile:
 *
 * Read a list of files from a file and add them to the filelist for
 * input data.  The filename is expected to be the last
 * space-separated field on the line.
 *
 * Returns the number of file names parsed from the list or -1 on error.
 ***************************************************************************/
static int
readlistfile (char *listfile)
{
  FILE *fp;
  char  line[1024];
  char *ptr;
  int   filecnt = 0;
  
  char  filename[1024];
  char *lastfield = 0;
  int   fields = 0;
  int   wspace;
  
  /* Open the list file */
  if ( (fp = fopen (listfile, "rb")) == NULL )
    {
      if (errno == ENOENT)
        {
          fprintf (stderr, "Could not find list file %s\n", listfile);
          return -1;
        }
      else
        {
          fprintf (stderr, "Error opening list file %s: %s\n",
		   listfile, strerror (errno));
          return -1;
        }
    }
  
  if ( verbose )
    fprintf (stderr, "Reading list of input files from %s\n", listfile);
  
  while ( (fgets (line, sizeof(line), fp)) !=  NULL)
    {
      /* Truncate line at first \r or \n, count space-separated fields
       * and track last field */
      fields = 0;
      wspace = 0;
      ptr = line;
      while ( *ptr )
	{
	  if ( *ptr == '\r' || *ptr == '\n' || *ptr == '\0' )
	    {
	      *ptr = '\0';
	      break;
	    }
	  else if ( *ptr != ' ' )
	    {
	      if ( wspace || ptr == line )
		{
		  fields++; lastfield = ptr;
		}
	      wspace = 0;
	    }
	  else
	    {
	      wspace = 1;
	    }
	  
	  ptr++;
	}
      
      /* Skip empty lines */
      if ( ! lastfield )
	continue;
      
      if ( fields >= 1 && fields <= 3 )
	{
	  fields = sscanf (lastfield, "%s", filename);
	  
	  if ( fields != 1 )
	    {
	      fprintf (stderr, "Error parsing file name from: %s\n", line);
	      continue;
	    }
	  
	  if ( verbose > 1 )
	    fprintf (stderr, "Adding '%s' to input file list\n", filename);
	  
	  addnode (&filelist, NULL, filename);
	  filecnt++;
	  
	  continue;
	}
    }
  
  fclose (fp);
  
  return filecnt;
}  /* End readlistfile() */


/***************************************************************************
 * addnode:
 *
 * Add node to the specified list.
 ***************************************************************************/
static void
addnode (struct listnode **listroot, char *key, char *data)
{
  struct listnode *lastlp, *newlp;
  
  if ( data == NULL )
    {
      fprintf (stderr, "addnode(): No file name specified\n");
      return;
    }
  
  lastlp = *listroot;
  while ( lastlp != 0 )
    {
      if ( lastlp->next == 0 )
        break;
      
      lastlp = lastlp->next;
    }
  
  newlp = (struct listnode *) malloc (sizeof (struct listnode));
  memset (newlp, 0, sizeof (struct listnode));
  if ( key ) newlp->key = strdup(key);
  else newlp->key = key;
  if ( data) newlp->data = strdup(data);
  else newlp->data = data;
  newlp->next = 0;
  
  if ( lastlp == 0 )
    *listroot = newlp;
  else
    lastlp->next = newlp;
  
}  /* End of addnode() */


/***************************************************************************
 * record_handler:
 * Saves passed records to the output file.
 ***************************************************************************/
static void
record_handler (char *record, int reclen, void *vofp)
{
  if ( fwrite(record, reclen, 1, (FILE *)vofp) != 1 )
    {
      fprintf (stderr, "Error writing to output file\n");
    }
}  /* End of record_handler() */


/***************************************************************************
 * usage:
 * Print the usage message and exit.
 ***************************************************************************/
static void
usage (void)
{
  fprintf (stderr, "%s version: %s\n\n", PACKAGE, VERSION);
  fprintf (stderr, "Convert MT bin file time series data to Mini-SEED.\n\n");
  fprintf (stderr, "Usage: %s [options] file1 [file2 file3 ...]\n\n", PACKAGE);
  fprintf (stderr,
	   " ## Options ##\n"
	   " -V             Report program version\n"
	   " -h             Show this usage message\n"
	   " -v             Be more verbose, multiple flags can be used\n"
	   " -S             Include SEED blockette 100 for very irrational sample rates\n"
	   " -C             Create a separate output file for each channel segment\n"
	   " -n network     Specify the SEED network code (currently %s)\n"
	   " -s station     Specify the SEED station code, default is blank\n"
	   " -l location    Specify the SEED location code, default is blank\n"
	   " -r bytes       Specify record length in bytes for packing, default: 4096\n"
	   " -e encoding    Specify SEED encoding format for packing, default: 11 (Steim2)\n"
	   " -b byteorder   Specify byte order for packing, MSBF: 1 (default), LSBF: 0\n"
	   "\n"
	   " -o outfile     Specify output file, default is %s.STA.yyyy-mm-ddTHH:MM:SS\n"
	   "\n"
	   " file(s)        File(s) of input data\n"
	   "                  If a file is prefixed with an '@' it is assumed to contain\n"
	   "                  a list of data files to be read\n"
	   "\n"
	   "Supported Mini-SEED encoding formats:\n"
           " 3  : 32-bit integers\n"
           " 10 : Steim 1 compression 32-bit integers\n"
           " 11 : Steim 2 compression 32-bit integers (default)\n"
	   "\n", network, network);
}  /* End of usage() */
