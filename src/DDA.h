/*
   DDA.h --
   	This header file defines global constants for the DDA program.
 */

#ifndef _DDA_H_
#define _DDA_H_

#define DDA_VERSION "0.9.0"

/* Length of an input token, e.g. "x:", "#", ... */
#define TOK_LEN 64
#define TOK_LEN_FMT "%63s"

/* Maximum length of a file or directory name */
#define PATH_LEN_MAX	4096
#define PATH_LEN_FMT	"%4095s"

/* Maximum number of radars */
#define MAX_RADARS 3

/* Global functions. See definitions for more information */
int DDA_Get_Radar_Idx(const char *);
char *DDA_Get_Radar_Nm(int);

#endif
