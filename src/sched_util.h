#ifndef SCHED_UTIL
#define SCHED_UTIL


struct abstract_element
{
  struct abstract_element *next;
  double t;
};


struct schedule_helper
{
  double dt;				/* timestep per slot */
  double dt_1;				/* dt_1 = 1/dt */
  double now;				/* start time of the scheduler */
  int buf_len;				/* number of slots in the scheduler */
  int index;				/* points to the next time block */
  int count;				/* number of future items */
  int *circ_buf_count;			/* array of counts of items that will
					   happen at time equal to the index 
					   value of the array */
  struct abstract_element **circ_buf_head;   /* array of the heads of the linked
						lists for future items */
  struct abstract_element **circ_buf_tail;   /* array of the tails of the linked
						lists for future items */
  struct abstract_element *current;	   /* points to the current item */
  struct abstract_element *current_tail;   /* points to the tail of the linked
					      list of the current items */
  int current_count;			   /* number of current items */
  struct schedule_helper *next_scale;
  int error;				   /* error code (1 - on error, 0 -
						no errors) */
};

struct abstract_element* ae_list_sort(struct abstract_element *ae);

struct schedule_helper* create_scheduler(double dt_min,double dt_max,int maxlen,double start_time);

int schedule_insert(struct schedule_helper *sh,void *data,int put_neg_in_current);
void schedule_excert(struct schedule_helper *sh,void *data,void *blank,int size);
int schedule_advance(struct schedule_helper *sh,struct abstract_element **head,
                     struct abstract_element **tail);
                     
void schedule_sort(struct schedule_helper *sh);
void* schedule_next(struct schedule_helper *sh);
#define schedule_add(x,y) schedule_insert((x),(y),1)

int schedule_anticipate(struct schedule_helper *sh,double *t);

void delete_scheduler(struct schedule_helper *sh);


#endif

