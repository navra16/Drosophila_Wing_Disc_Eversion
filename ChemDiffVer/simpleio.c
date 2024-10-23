/*
*				simplio:
*
*	Contains routines for formatting vectors, arrays etc:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include "cdecs.h"
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/utsname.h>

#if defined(float)
	static const char *scan_float_fmt = "%lf";
#else /* defined(float) */
	static const char *scan_float_fmt = "%f";
#endif /* defined(float) */


/*
*			print_title():
*
*	Prints a Run-time Header.
*/

EXPORT void print_title(
	FILE       *file,
	const char *title)
{
	time_t		tvec;
	
	(void) time(&tvec);		/* Get Run Time */

	(void) fprintf(file,"%s",title);
	(void) fprintf(file,"\n\t\tDATE OF RUN  %s\n\n",ctime(&tvec));
	print_machine_parameters(file);
}		/*end print_title*/

EXPORT	void	print_machine_parameters(
	FILE *file)
{
	struct utsname Uts;

	(void) foutput(file);
	(void) fprintf(file,"#MACHINE PARAMETERS\n");
	(void) uname(&Uts);
	(void) fprintf(file,"#\tHostname                 = %s\n",Uts.nodename);
	(void) fprintf(file,"#\tOperating System         = %s\n",Uts.sysname);
	(void) fprintf(file,"#\tOS Release               = %s\n",Uts.release);
	(void) fprintf(file,"#\tOS Version               = %s\n",Uts.version);
	(void) fprintf(file,"#\tCPU Type                 = %s\n",Uts.machine);
	(void) fprintf(file,"#\tByte Ordering            = ");
	switch (ft_endian_type())
	{
	case FT_BIG_ENDIAN:
	    (void) fprintf(file,"Big Endian\n");
	    break;
	case FT_LITTLE_ENDIAN:
	    (void) fprintf(file,"Little Endian\n");
	    break;
	case FT_UNKNOWN_ENDIAN:
	default:
 	    (void) printf("Undetermined Endian\n");
	    break;
	}
	(void) fprintf(file,"#\tFloating Point Word Size = %lu\n",
	               sizeof(float));
	(void) fprintf(file,"\n");
}		/*end print_machine_parameters*/

EXPORT	void	stripcomm(
	char	   *scfname,	/*Strip commented file name*/
	const char *fname)	/*Raw input file name*/
{
	FILE	*file = fopen(fname,"r");
	FILE	*scfile;
	char	*c, line[2048];
	static const char *separator = ": ";
	bool	status;
	size_t	sep_len = strlen(separator);
	long	io_pid;

	if (file == NULL)
	{
	    screen("ERROR in stripcomm(), can't open %s\n",fname);
	    clean_up(ERROR);
	}
	if (fgetstring(file,separator) == FUNCTION_FAILED)
	{
	    /*File is already stripcommented*/
	    (void) strcpy(scfname,fname);
	    (void) fclose(file);
	    return;
	}

	io_pid = (is_io_node(pp_mynode())) ? (long) getpid() : 0;
	pp_global_lmax(&io_pid,1L);

	(void) sprintf(scfname,"%s-%ld.sc",fname,io_pid);

	if (is_io_node(pp_mynode()))
	{
	    rewind(file);
	    if ((scfile = fopen(scfname,"w")) == NULL)
	    	status = NO;
	    else
	    {
	    	status = YES;
	    	while (fgets(line,2046,file) != NULL)
	    	    if ((c = strstr(line,separator)) != NULL)
		    	(void) fprintf(scfile,"%s",c+sep_len);
	        (void) fclose(scfile);
	    }
	}
	else
	    status = YES;

	(void) fclose(file);

	if (pp_min_status(status) == NO)
	{
	    screen("ERROR in stripcomm(), can't open %s\n",scfname);
	    clean_up(ERROR);
	}

	return;
}		/*end stripcomm*/

EXPORT	uint64_t u_ptr2ull(void* p)
{
	uint64_t ll = (unsigned long)p;
	return ll;
}		/*end u_ptr2ull*/

/*
*				fprint_vector():
*
*	Prints an n-component vector with a title at the top:
*	The vector is a vector of floating point numbers.
*	The numbers are printed according to the given format
*	with no blank between them them.   The format may
*	include blanks or newlines.
*	Usage:
*			fprint_vector(file,title,n,vector,format)
*			FILE *file;
*			const char *title;
*			int n;
*			float *vector;
*			const char *format;
*/


EXPORT void fprint_vector(
	FILE	   *file,
	const char *title,
	int	   n,
	float	   *vector,
	const char *format)
{
	int		i;


	(void) fprintf(file,"\n\n");
	if (title != NULL) (void) fprintf(file,"%s\n",title);
	for (i = 0; i < n; i++)
	    (void) fprintf(file,format,vector[i]);
	(void) fprintf(file,"\n\n");
}		/*end fprint_vector*/









/*
*				fprint_vector_as_matrix():
*
*	Prints a float vector of length len as a matrix of 'cols' 
*	columns.   Individual numbers are printed according to the
*	given format with no blank between them.   The format may
*	include extra blanks.   The given title is printed at the
*	top.   
*/

EXPORT void fprint_vector_as_matrix(
	FILE		*file,
	const char	*title,
	int		length,
	int		cols,
	float		*vector,
	const char	*format)
{
	int		i,col;
	bool		KeepPrinting = YES;

	(void) fprintf(file,"\n\n");
	if (title!=NULL) (void) fprintf(file,"%s\n",title);
	i = 0;
	while (KeepPrinting == YES)
	{
	    for (col = 0; i < length && col < cols; col++,i++)
	    	(void) fprintf(file,format,vector[i]);
	    (void) fprintf(file,"\n");
	    if (i == length)
		break;
	}
	(void) fprintf(file,"\n\n");
}		/*end fprint_vector_as_matrix*/









/*
*				fprint_matrix():
*
*	Prints a float matrix of 'rows' rows and 'cols'  columns.
*	Individual numbersin the matrix are printed according to the
*	given format with no blank between them.   The format may
*	include extra blanks.   The given title is printed at the
*	top.
*/

EXPORT  void    fprint_matrix(
	FILE		*file,
	const char	*title,
	int		rows,
	int		cols,
	float		**matrix,
	const char      *format)
{
	int		row,col;

	(void) fprintf(file,"\n\n");
	if (title != NULL) (void) fprintf(file,"%s\n",title);
	for (row = rows-1; row >= 0; row--)
	{
	    for(col=0; col<cols; col++)
	    	(void) fprintf(file,format,matrix[row][col]);
	    (void) fprintf(file,"\n");
	}
	(void) fprintf(file,"\n\n");
}		/*end fprint_matrix*/





/*
*				fprint_float_vector():
*
*	Prints an n-component vector with a title at the top:
*	The vector is a vector of floating point numbers.
*	The numbers are printed according to the given format
*	with no blank between them them.   The format may
*	include blanks or newlines.
*	Usage:
*			fprint_float_vector(file,title,n,vector,format)
*			FILE *file;
*			const char *title;
*			int n;
*			float *vector;
*			const char *format;
*/


EXPORT  void    fprint_float_vector(
	FILE		*file,
	const char	*title,
	int		n,
	float		*vector,
	const char	*format)
{
	int		i;


	(void) fprintf(file,"\n\n");
	if (title != NULL) (void) fprintf(file,"%s\n",title);
	for (i = 0; i < n; i++)
	    (void) fprintf(file,format,vector[i]);
	(void) fprintf(file,"\n\n");
}		/*end fprint_float_vector*/









/*
*			fprint_float_vector_as_matrix():
*
*	Prints a float vector of length len as a matrix of 'cols' 
*	columns.   Individual numbers are printed according to the
*	given format with no blank between them.   The format may
*	include extra blanks.   The given title is printed at the
*	top.   
*/

EXPORT  void    fprint_float_vector_as_matrix(
	FILE		*file,
	const char      *title,
	int		length,
	int		cols,
	float		*vector,
	const char      *format)
{
	int		i,col;
	bool		KeepPrinting = YES;

	(void) fprintf(file,"\n\n");
	if (title != NULL) (void) fprintf(file,"%s\n",title);
	i = 0;
	while (KeepPrinting == YES)
	{
	    for (col = 0; i<length && col < cols; col++,i++)
	    	(void) fprintf(file,format,vector[i]);
	    (void) fprintf(file,"\n");
	    if (i == length)
		break;
	}
	(void) fprintf(file,"\n\n");
}		/*end fprint_float_vector_as_matrix*/









/*
*			fprint_float_matrix():
*
*	Prints a float matrix of 'rows' rows and 'cols'  columns.
*	Individual numbers in the matrix are printed according to the
*	given format with no blank between them.   The format may
*	include extra blanks.   The given title is printed at the
*	top.
*/

EXPORT  void    fprint_float_matrix(
	FILE		*file,
	const char      *title,
	int		rows,
	int		cols,
	float		**matrix,
	const char      *format)
{
	int		row,col;

	(void) fprintf(file,"\n\n");
	if (title != NULL) (void) fprintf(file,"%s\n",title);
	for (row = rows-1; row >= 0; row--)
	{
	    for(col=0; col<cols; col++)
	    	(void) fprintf(file,format,matrix[row][col]);
	    (void) fprintf(file,"\n");
	}
	(void) fprintf(file,"\n\n");
}		/*end fprint_float_matrix*/







/*
*			fprint_double_vector():
*
*	Prints an n-component vector with a title at the top:
*	The vector is a vector of double precision point numbers.
*	The numbers are printed according to the given format
*	with no blank between them them.   The format may
*	include blanks or newlines.
*	Usage:
*			fprint_double_vector(file,title,n,vector,format)
*			FILE *file;
*			const char *title;
*			int n;
*			double *vector;
*			const chart *format;
*/


EXPORT  void    fprint_double_vector(
	FILE		*file,
	const char      *title,
	int		n,
	double		*vector,
	const char      *format)
{
	int		i;


	(void) fprintf(file,"\n\n");
	if (title != NULL)
	    (void) fprintf(file,"%s\n",title);
	for (i = 0; i < n; i++)
	    (void) fprintf(file,format,vector[i]);
	(void) fprintf(file,"\n\n");
}		/*end fprint_double_vector*/









/*
*			fprint_double_vector_as_matrix():
*
*	Prints a double vector of length len as a matrix of 'cols' 
*	columns.   Individual numbers are printed according to the
*	given format with no blank between them.   The format may
*	include extra blanks.   The given title is printed at the
*	top.   
*/

EXPORT  void    fprint_double_vector_as_matrix(
	FILE		*file,
	const char      *title,
	int		length,
	int		cols,
	double		*vector,
	const char      *format)
{
	int		i,col;
	bool		KeepPrinting = YES;

	(void) fprintf(file,"\n\n");
	if (title != NULL)
	    (void) fprintf(file,"%s\n",title);
	i = 0;
	while (KeepPrinting == YES)
	{
	    for (col=0; i<length && col < cols; col++,i++)
	    	(void) fprintf(file,format,vector[i]);
	    (void) fprintf(file,"\n");
	    if (i == length)
		break;
	}
	(void) fprintf(file,"\n\n");
}		/*end fprint_double_vector_as_matrix*/









/*
*			fprint_double_matrix():
*
*	Prints a double matrix of 'rows' rows and 'cols'  columns.
*	Individual numbersin the matrix are printed according to the
*	given format with no blank between them.   The format may
*	include extra blanks.   The given title is printed at the
*	top.
*/

EXPORT  void    fprint_double_matrix(
	FILE	   *file,
	const char *title,
	int	   rows,
	int	   cols,
	double	   **matrix,
	const char *format)
{
	int		row,col;

	(void) fprintf(file,"\n\n");
	if (title != NULL) (void) fprintf(file,"%s\n",title);
	for (row = rows-1; row >= 0; row--)
	{
	    for (col = 0; col < cols; col++)
	    	(void) fprintf(file,format,matrix[row][col]);
	    (void) fprintf(file,"\n");
	}
	(void) fprintf(file,"\n\n");
}		/*end fprint_double_matrix*/




/*
*			fprint_int_vector_as_matrix():
*
*	Prints a integer vector of length len as a matrix of 'cols' 
*	columns.   Individual numbers are printed according to the
*	given format with no blank between them.   The format may
*	include extra blanks.   The given title is printed at the
*	top.   
*/

EXPORT  void    fprint_int_vector_as_matrix(
	FILE		*file,
	const char      *title,
	int		length,
	int		cols,
	int		*vector,
	const char      *format)
{
	int		i,col;
	bool		KeepPrinting = YES;

	(void) fprintf(file,"\n\n");
	if (title != NULL) (void) fprintf(file,"%s\n",title);
	i = 0;
	while (KeepPrinting == YES)
	{
	    for (col = 0; i < length && col < cols; col++,i++)
	    	(void) fprintf(file,format,vector[i]);
	    (void) fprintf(file,"\n");
	    if (i == length)
		break;
	}
	(void) fprintf(file,"\n\n");
}		/*end fprint_int_vector_as_matrix*/






/*
*			fprint_int_matrix():
*
*	Prints a integer matrix of 'rows' rows and 'cols'  columns.
*	Individual numbers in the matrix are printed according to the
*	given format with no blank between them.   The format may
*	include extra blanks.   The given title is printed at the
*	top.
*/

EXPORT  void    fprint_int_matrix(
	FILE		*file,
	const char      *title,
	int		rows,
	int		cols,
	int		**matrix,
	const char      *format)
{
	int		row,col;

	(void) fprintf(file,"\n\n");
	if (title != NULL) (void) fprintf(file,"%s\n",title);
	for (row = rows-1; row >= 0; row--)
	{
	    for (col = 0; col < cols; col++)
	    	(void) fprintf(file,format,matrix[row][col]);
	    (void) fprintf(file,"\n");
	}
	(void) fprintf(file,"\n\n");
}		/*end fprint_int_matrix*/



LOCAL bool binary_output = YES;

EXPORT void set_binary_output(
	bool		yes_no)
{
	binary_output = yes_no;
}		/*end set_binary_output*/

EXPORT bool is_binary_output(void)
{
	return binary_output;
}		/*end is_binary_output*/

EXPORT	void	fwrite_float(
	FILE	   *file,
	const char *msg,
	float	   x,
	bool	   bio,
	const char *fmt,
	const char *end)
{
	if (msg != NULL)
	    (void) fprintf(file,"%s",msg);
	if (bio == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) &x,sizeof(float),1,file);
	}
	else
	    (void) fprintf(file,fmt,x);
	(void) fprintf(file,"%s",end);
}		/* end fwrite_float*/

EXPORT	size_t	read_binary_real_array(
	float         *a,
	size_t        len,
	const IO_TYPE *io_type)
{
	size_t nread = 0;
	FILE *file = io_type->file;
	size_t  i;

	if (io_type->read_float_size == io_type->cpu_float_size)
	{
	    nread = fread((void *) a,sizeof(float),len,file);
	    if (io_type->reverse_endian)
	        for (i = 0; i < nread; ++i)
		    byte_reverse(a[i]);
	}
	else if (io_type->read_float_size > io_type->cpu_float_size)
	{
	    double tmp;
	    for (i = 0; i < len; ++i)
	    {
	    	if (fread((void *) &tmp,sizeof(double),1,file))
		{
		    ++nread;
		    if (io_type->reverse_endian)
		        byte_reverse(tmp);
	    	    a[i] = tmp;
	        }
	    }
	}
	else
	{
	    truefloat tmp;
	    for (i = 0; i < len; ++i)
	    {
	    	if (fread((void *) &tmp,sizeof(truefloat),1,file))
		{
		    ++nread;
		    if (io_type->reverse_endian)
		        byte_reverse(tmp);
	    	    a[i] = tmp;
	        }
	    }
	}
	return nread;
}		/*end read_binary_real_array*/

EXPORT	size_t	read_binary_int_array(
	int           *a,
	size_t        len,
	const IO_TYPE *io_type)
{
	FILE   *file = io_type->file;
	size_t i;
	size_t nread = 0;

	nread = fread((void *) a,sizeof(int),len,file);
	if (io_type->reverse_endian)
	    for (i = 0; i < nread; ++i)
		byte_reverse(a[i]);
	return nread;
}		/*end read_binary_int_array*/

EXPORT	size_t	read_binary_uint_array(
	unsigned int  *a,
	size_t        len,
	const IO_TYPE *io_type)
{
	FILE   *file = io_type->file;
	size_t i;
	size_t nread = 0;

	nread = fread((void *) a,sizeof(unsigned int),len,file);
	if (io_type->reverse_endian)
	    for (i = 0; i < nread; ++i)
		byte_reverse(a[i]);
	return nread;
}		/*end read_binary_uint_array*/

EXPORT	size_t	read_binary_uint64_t_array(
	uint64_t      *a,
	size_t        len,
	const IO_TYPE *io_type)
{
	FILE   *file = io_type->file;
	size_t i;
	size_t nread = 0;

	nread = fread((void *) a,sizeof(uint64_t),len,file);
	if (io_type->reverse_endian)
	    for (i = 0; i < nread; ++i)
		byte_reverse(a[i]);
	return nread;
}		/*end read_binary_uint_array*/


EXPORT	float	fread_float(
	const char    *label,
	const IO_TYPE *io_type)
{
	FILE	*file = io_type->file;
	float	x;
	char	s[81];
	int	c;

	if (label != NULL)
	    (void) fgetstring(file,label);
	if ((c = getc(file)) != '\f') /*NOBINARY*/
	{
	    long current;
	    (void) ungetc(c,file);
	    current = ftell(file);
	    (void) fgets(s,5,file);
	    if (strncasecmp("inf",s,3) == 0)
	    	x = HUGE_VAL;
	    else if (strncasecmp("-inf",s,4) == 0)
	    	x = -HUGE_VAL;
	    else
	    {
		(void) fseek(file,current,SEEK_SET);
	        if (fscanf(file,scan_float_fmt,&x) != 1)
		{
	    	    screen("ERROR in fread_float(), "
	    	           "can't read floating point variable\n");
	    	    clean_up(ERROR);
		}
	    }
	}
	else
	{
	    (void) getc(file);
	    if (io_type->read_float_size == io_type->cpu_float_size)
	    {
	        (void) fread((void *)&x,sizeof(float),1,file);
	        if (io_type->reverse_endian)
	            byte_reverse(x);
	    }
	    else if (io_type->read_float_size > io_type->cpu_float_size)
	    {
		double tmp;
	        (void) fread((void *)&tmp,sizeof(double),1,file);
	        if (io_type->reverse_endian)
	            byte_reverse(tmp);
		x = tmp;
	    }
	    else
	    {
	        truefloat tmp;
	        (void) fread((void *)&tmp,sizeof(truefloat),1,file);
	        if (io_type->reverse_endian)
	            byte_reverse(tmp);
		x = tmp;
	    }
	}
	return x;
}		/* end fread_float*/



/*
*			fprint_vector_of_floats():
*
*	Prints the first number elements of the vector f[] of floats
*	in binary or non-binary format according to the return value
*	of binary_output(). 
*
*	In binary mode, numbers are written out in the form:
*		\fnf[1]f[2]f[3].....
*	where n is the number of floats as a character.   In case
*	number > 255, the floats beyond 255 are printed using further
*	\fn sequence until all points are exhausted.
*
*	In non-binary mode numbers are written in a line, blank spaced
*	and with a newline at the end.
*/

EXPORT void fprint_vector_of_floats(
	FILE		*file,
	int		number,
	float		*f)
{
	int		i, num;
	static const	int	MAX_FLOATS = 255;

	if (is_binary_output() == YES)
	{
	    while (number > 0)
	    {
	    	num = min(number,MAX_FLOATS);
	    	(void) fprintf(file,"\f%c",num);
	    	(void) fwrite((const void *)f,sizeof(float),num,file);
	    	number -= num;
	    }
	}
	else
	{
	    (void) fprintf(file,"%"FFMT,f[0]);
	    for (i = 1; i < number; i++)
	    	(void) fprintf(file," %"FFMT,f[i]);
	    (void) fprintf(file,"\n");
	}
}		/*end fprint_vector_of_floats*/

EXPORT	int	sscan_float(
	const char *s,
	float	   *px)
{
	int n;

	if (s[0] == '\0')
	    return 0;

	n = 1;
	if (strncasecmp("inf",s,3) == 0)
	    *px = HUGE_VAL;
	else if (strncasecmp("-inf",s,4) == 0)
	    *px = -HUGE_VAL;
	else
	    n = sscanf(s,scan_float_fmt,px);
	return n;
}		/*end sscan_float*/

EXPORT	int	fscan_float(
	FILE	*file,
	float	*px)
{
	char s[256];
	int  n;

	n = fscanf(file,"%s",s);
	if (n != 1)
	    return 0;
	if (strncasecmp("inf",s,3) == 0)
	    *px = HUGE_VAL;
	else if (strncasecmp("-inf",s,4) == 0)
	    *px = -HUGE_VAL;
	else
	    n = sscanf(s,scan_float_fmt,px);
	return n;
}		/*end fscan_float*/

EXPORT	bool	fread_bool(
	FILE	*file)
{
	bool    value;
	int     c;
	int	value_read;

	value = NO;
	c = getc(file);
	if (isalpha(c))
	{
	    char s[120];
	    (void) ungetc(c,file);
	    (void) fscanf(file,"%s",s);
	    if ((c == 'y') || (c == 'Y'))
	        value = YES;
	    else if ((c == 'n') || (c == 'N'))
	        value = NO;
	    else
	    {
	        screen("ERROR in fread_bool(), invalid value %s\n",s);
		clean_up(ERROR);
	    }
	}
	else
	{
	    (void) fscanf(file,"%d",&value_read);
	    switch (value_read)
	    {
	    case YES:
	        value = YES;
	    case NO:
	    default:
	        value = NO;
	    }
	}
	return value;
}		/*end fread_bool*/

EXPORT	void fprint_float(
	FILE       *file,
	const char *mesg,
	float      value,
	const char *end)
{
	(void) fprintf(file,"%s",mesg);
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) &value,sizeof(float),1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT,value);
	
	(void) fprintf(file,"%s",end);
}		/*end fprint_float*/

EXPORT	float	read_print_float(
	const char    *search_string,
	float	      dflt_value,
	const IO_TYPE *io_type)
{
	FILE	*file = io_type->file;
	float	value = dflt_value;

	if ((search_string != NULL) && (!fgetstring(file,search_string)))
	    return value; 
	value = fread_float(NULL,io_type);
	return value;
}		/*end read_print_float*/
