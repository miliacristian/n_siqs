#include "timing.h"
struct timespec timer;

void gettime(struct timespec*timer){
	if(timer==NULL){
		handle_error_with_exit("error in gettime\n");
	}
	if(clock_gettime(CLOCK_MONOTONIC,(struct timespec*)timer)!=0){
        	handle_error_with_exit("error in clock gettime\n");
    	}
}

struct timespec diff_timespec(struct timespec time_current,struct timespec timer){
	//timecurrent>timer
	struct timespec time_sub;
	if(time_current.tv_nsec>=timer.tv_nsec){
		time_sub.tv_nsec=time_current.tv_nsec-timer.tv_nsec;
		time_sub.tv_sec=time_current.tv_sec-timer.tv_sec;
	}
	else{//riporto nanosecondi di time_current minori di nanosecondi time
		time_current.tv_sec-=1;
		timer.tv_nsec-=time_current.tv_nsec;
		time_current.tv_nsec=0;//evita l'overflow nei nanosecondi
		time_current.tv_nsec+=1000000000;
		time_sub.tv_nsec=time_current.tv_nsec-timer.tv_nsec;
		time_sub.tv_sec=time_current.tv_sec-timer.tv_sec;
	}
	return time_sub;
}
void print_time_elapsed_local(char*string,struct timespec*timer_thread,unsigned int tid){
	//ns=1000000000=1 sec
	//ns=1000000 1 ms
	//calculate timer relative to timer_thread,change variable timer_thread to actual timer taken by gettime
	long ns,ms,sec,min,hour,temp;
	if(string==NULL || timer_thread==NULL){
		handle_error_with_exit("error in print time elapsed local\n");
	}
	struct timespec time_current,time_sub;
	if(clock_gettime(CLOCK_MONOTONIC,&time_current)!=0){
		handle_error_with_exit("error in print_time\n");
	}
	time_sub=diff_timespec(time_current,*timer_thread);
	if(clock_gettime(CLOCK_MONOTONIC,timer_thread)!=0){
		handle_error_with_exit("error in clockgettime\n");
	}
	ns=time_sub.tv_nsec%1000000;
	time_sub.tv_nsec-=ns;
	ms=time_sub.tv_nsec/1000000;
	sec=time_sub.tv_sec%60;//i secondi sono modulo 60
	min=(time_sub.tv_sec-sec)/60;
	temp=min%60;//temp min
	hour=(min-temp)/60;
	min=temp;
	printf("%s tid %u:hour=%ld min=%ld sec:%ld ms=%ld ns:%ld\n",string,tid,hour,min,sec,ms,ns);
	return;
}
void print_total_time_elapsed(char*string,struct timespec timer_start){
	//calculate timer relative to timer_start,don't change any variable
	long ns,ms,sec,min,hour,temp;
	if(string==NULL){
		handle_error_with_exit("error in print time\n");
	}
	struct timespec time_temp,time_sub;
	if(clock_gettime(CLOCK_MONOTONIC,&time_temp)!=0){
        	handle_error_with_exit("error in print_time\n");
    	}
	time_sub=diff_timespec(time_temp,timer_start);
	ns=time_sub.tv_nsec%1000000;
	time_sub.tv_nsec-=ns;
	ms=time_sub.tv_nsec/1000000;
	sec=time_sub.tv_sec%60;//i secondi sono modulo 60
	min=(time_sub.tv_sec-sec)/60;
	temp=min%60;//temp min
	hour=(min-temp)/60;
	min=temp;
	printf("%s:hour=%ld min=%ld sec:%ld ms=%ld ns:%ld\n",string,hour,min,sec,ms,ns);
	return;
}

struct timespec print_time_elapsed(char*string){
	//ns=1000000000=1 sec
	//ns=1000000 1 ms
	//change global timer variable,this function is not thread_safe
	long ns,ms,sec,min,hour,temp;
	if(string==NULL){
		handle_error_with_exit("error in print time\n");
	}
	struct timespec time_current,time_sub;
	if(clock_gettime(CLOCK_MONOTONIC,&time_current)!=0){
        	handle_error_with_exit("error in print_time\n");
    	}
	time_sub=diff_timespec(time_current,timer);
	//re init timer
	timer=time_current;
	ns=time_sub.tv_nsec%1000000;
	time_sub.tv_nsec-=ns;
	ms=time_sub.tv_nsec/1000000;
	sec=time_sub.tv_sec%60;//i secondi sono modulo 60
	min=(time_sub.tv_sec-sec)/60;
	temp=min%60;//temp min
	hour=(min-temp)/60;
	min=temp;
	printf("%s:hour=%ld min=%ld sec:%ld ms=%ld ns:%ld\n",string,hour,min,sec,ms,ns);
	return time_sub;
}
