/*	$Id: sw-vector.h,v 1.7 2009/06/16 23:26:21 rumble Exp $	*/

int sw_vector_cleanup(void);
int	sw_vector_setup(int, int, int, int, int, int, int, int, int, bool);
void	sw_vector_stats(uint64_t *, uint64_t *, double *);
int	sw_vector(char const *, int, int, char const *, int, char const *, int, bool);
