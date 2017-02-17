#ifndef THREADING_H_
#define THREADING_H_

#include <iostream>
#include "spinlock.h"

// Note that USE_SPINLOCK trumps BOWTIE_PTHREADS

#ifdef BOWTIE_PTHREADS
#include <pthread.h>
#endif

#ifdef USE_SPINLOCK
#  include "spinlock.h"
#  define MUTEX_T SpinLock
#  define MUTEX_INIT(l)
#  define MUTEX_LOCK(l) (l).Enter()
#  define MUTEX_UNLOCK(l) (l).Leave()
#else
#  ifdef BOWTIE_PTHREADS
#    define MUTEX_T pthread_mutex_t
#    define MUTEX_INIT(l) pthread_mutex_init(&l, NULL)
#    define MUTEX_LOCK(l) pthread_mutex_lock(&l)
#    define MUTEX_UNLOCK(l) pthread_mutex_unlock(&l)
#  else
#    define MUTEX_T int
#    define MUTEX_INIT(l) l = 0
#    define MUTEX_LOCK(l) l = 1
#    define MUTEX_UNLOCK(l) l = 0
#  endif /* BOWTIE_PTHREADS */
#endif /* USE_SPINLOCK */

#ifdef BOWTIE_PTHREADS
static inline void joinThread(pthread_t th) {
	int ret, *tmp;
	if((ret = pthread_join(th, (void**)(int**)&tmp)) != 0) {
		std::cerr << "Error: pthread_join returned non-zero status: "
		          << ret << std::endl;
		throw 1;
	}
}

static inline void createThread(pthread_t* th,
                                void *(*start_routine) (void *),
                                void *arg)
{
	int ret;
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	pthread_attr_setstacksize(&pt_attr, 2 << 20);
	if((ret = pthread_create(th, &pt_attr, start_routine, arg)) != 0) {
		std::cerr << "Error: pthread_create returned non-zero status: "
		          << ret << std::endl;
		throw 1;
	}
}
#endif

/**
 * Wrap a lock; obtain lock upon construction, release upon destruction.
 */
class ThreadSafe {
public:
	ThreadSafe(MUTEX_T* lock) {
		lock_ = lock;
		MUTEX_LOCK(*lock_);
	}
	~ThreadSafe() { MUTEX_UNLOCK(*lock_); }
private:
	MUTEX_T *lock_;
};

#endif
