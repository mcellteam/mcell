#ifndef SCHED_UTIL
#define SCHED_UTIL


struct abstract_element
{
  struct abstract_element *next;
  double t;
};


struct schedule_helper
{
  double dt,dt_1,now;
  int buf_len,index,count;
  int *circ_buf_count;
  struct abstract_element **circ_buf_head;
  struct abstract_element **circ_buf_tail;
  struct abstract_element *current;
  struct abstract_element *current_tail;
  int current_count;
  struct schedule_helper *next_scale;
};

struct abstract_element* ae_list_sort(struct abstract_element *ae);

struct schedule_helper* create_scheduler(double dt_min,double dt_max,int maxlen,double start_time);

void schedule_insert(struct schedule_helper *sh,void *data,int put_neg_in_current);
void schedule_excert(struct schedule_helper *sh,void *data,void *blank,int size);
int schedule_advance(struct schedule_helper *sh, void** head, void** tail);

void schedule_sort(struct schedule_helper *sh);
void* schedule_next(struct schedule_helper *sh);
#define schedule_add(x,y) schedule_insert((x),(y),1)

int schedule_anticipate(struct schedule_helper *sh,double *t);

void delete_scheduler(struct schedule_helper *sh);


#endif

