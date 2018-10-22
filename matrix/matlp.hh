#ifndef MATLP__
#define MATLP__

#include <new>
#include "ross.h"

typedef tw_stime Time;

static const int nmat = 5;
static const double mean_delay = 1.0;
static const int debug = 0;

template<typename myuint>
myuint intinv(myuint b) {
  myuint t0 = 0,t = 1,q,r;
  myuint a = 0,nota = ~a;

  if(b <= 1) return b;

  q = nota / b;
  if(b*q+b == 0) return 0;
  r = a - q*b;
  while(r > 0) {
    const myuint temp = t0 - q*t;
    t0 = t;
    t = temp;
    a = b;
    b = r;
    q = a/b;
    r = a - q*b;
  }
  if(b == 1) return t;
  else return 0;
}

template<typename myuint>
int matinv(int n,myuint X[],myuint Y[]) {
  int p[n],q[n];

  for(int i = 0; i<n; i++)
    p[i] = q[i] = i;
    
  for(int i = 0; i<n; i++) {
    for(int j = 0; j<n; j++)
      Y[i*n+j] = 0;
    Y[i*n+i] = 1;
  }

  for(int i = 0; i<n; i++) {
    myuint xi = 0;
    /* Find pivot element */ {
      int j,k;
      for(j = i; j<n; j++)
	for(k = i; k<n; k++) {
	  xi = intinv(X[p[j]*n+q[k]]);
	  if(xi != 0) {
	    int t = p[i];
	    p[i] = p[j];
	    p[j] = t;
	    t = q[i];
	    q[i] = q[k];
	    q[k] = t;
	    goto found;
	  }
	}
    found:
      if(xi == 0) return 0;
    }
    
    /* Make diagonal element 1 */ {
      for(int j = i; j<n; j++)
	X[p[i]*n+q[j]] *= xi;
      for(int j = 0; j<n; j++)
	Y[p[i]*n+q[j]] *= xi;
    }

    /* Elimination */ {
      for(int j = i+1; j<n; j++) {
	myuint xji = X[p[j]*n+q[i]];
	for(int k = i; k<n; k++)
	  X[p[j]*n+q[k]] -= xji*X[p[i]*n+q[k]];
	for(int k = 0; k<n; k++)
	  Y[p[j]*n+k] -= xji*Y[p[i]*n+k];
      }
    }
  }

  /* Back substitution */ {
    for(int i = n-1; i>=0; i--) {
      for(int j = 0; j<i; j++) {
	myuint xji = X[p[j]*n+q[i]];
	for(int k = 0; k<n; k++)
	  Y[p[j]*n+k] -= xji*Y[p[i]*n+k];
	X[p[j]*n+q[i]] -= xji*X[p[i]*n+q[i]];
      }
    }
  }

  /* Undo permutation, and copy result to Y */ {
    for(int i = 0; i<n; i++)
      for(int j = 0; j<n; j++)
	X[q[i]*n+j] = Y[p[i]*n+j];
    for(int i = 0; i<n*n; i++)
      Y[i] = X[i];
  }
  return 1;
}

template<typename myuint>
void matmul(int n,myuint A[],myuint B[],myuint AB[]) {
  for(int i = 0; i<n; i++)
    for(int j = 0; j<n; j++) {
      myuint s = 0;
      for(int k = 0; k<n; k++)
	s = s + A[i*n+k]*B[k*n+j];
      AB[i*n+j] = s;
    }
}

template<typename myuint>
void matcopy(int n,myuint dest[],myuint src[]) {
  int n2 = n*n;
  for(int i = 0; i<n2; i++)
    dest[i] = src[i];
}

template<typename uint>
void matprint(int n,uint A[],const char *label) {
#if 1
  if(label)
    printf("  %s =\n",label);
  for(int i = 0; i<nmat; i++) {
    printf("    ");
    for(int j = 0; j<nmat; j++)
      printf("  %12llu",(unsigned long long int) A[i*nmat+j]);
    printf("\n");
  }
  printf("\n");
#endif
}

struct Event {
  typedef unsigned short uint;
  uint A[nmat*nmat];
  uint nextrand;
  int seq;
  void print(const char* label = 0) {
    matprint(nmat,A,label);
  }
};

class lp {
public:
  virtual void forward(tw_lp *twlp,Event *evt,const Time &t) = 0;
  virtual void backward(tw_lp *twlp,Event *evt,const Time &t) = 0;
  virtual void commit(tw_lp *twlp,Event *evt,const Time &t) = 0;
  virtual ~lp() {}

  static void forward_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp) {
    tw_stime tnow0 = tw_now(thislp);
    Time tnow = tnow0;
    lp* lpptr = static_cast<lp *>(state);
    Event* e = static_cast<Event *>(evt);
    lpptr->forward(thislp,e,tnow);
  }
  
  static void backward_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp) {
    Time tnow = tw_now(thislp);
    lp* lpptr = static_cast<lp *>(state);
    Event *e = static_cast<Event *>(evt);
    lpptr->backward(thislp,e,tnow);
  }
  
  static void commit_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp) {
    Time tnow = tw_now(thislp);
    lp* lpptr = static_cast<lp *>(state);
    Event* e = static_cast<Event *>(evt);
    lpptr->commit(thislp,e,tnow);
  }
};

template<typename myuint>
bool invcheck(int n,myuint A[]) {
  myuint Atmp[n*n],Ai[n*n];
  matcopy(n,Atmp,A);
  return matinv(n,Atmp,Ai) == 1;
}

template<typename myuint>
void rand_matrix(tw_lp *lpptr,int n,myuint A[],int imax) {
  do {
    for(int i = 0; i<n; i++)
      for(int j = 0; j<n; j++) {
	myuint x = tw_rand_integer(lpptr->rng,0,imax);
	A[i*n+j] = x;
      }
  } while(!invcheck(n,A));
}


struct matlp : public lp {
  typedef unsigned short uint;
  uint *M,*T; 
  const int n;

  struct global_stat_t {
    long long int nforward,nbackward,ncommit;
    uint A[nmat*nmat];
  };
  static global_stat_t global_stat;

  long long int nforward,nbackward,ncommit;
  int seq,lastseqcommit;


  void apply(uint X[]) {
    uint tmp[n*n];
    if(debug) {
      assert(invcheck(n,X));
      assert(invcheck(n,M));
    }

    matmul(n,M,X,tmp);
    matcopy(n,M,tmp);
  }

  void unapply(uint X[]) {
    uint Xtmp[n*n],Xi[n*n];

    if(debug)
      assert(invcheck(n,M));

    matcopy(n,Xtmp,X);
    assert(matinv(n,Xtmp,Xi) == 1);
    matmul(n,M,Xi,Xtmp);
    matcopy(n,M,Xtmp);
    
    if(debug)
      assert(invcheck(n,M));
  }


  matlp(tw_lp *lpptr,int nin = 2) : n(nin),nforward(0),nbackward(0),ncommit(0),seq(0),lastseqcommit(-1) {
    M = new uint[n*n];
    T = new uint[n*n];
    rand_matrix(lpptr,n,M,1000);
    rand_matrix(lpptr,n,T,1000);

    const int nevent_init = 10;
    for(int i = 0; i<nevent_init; i++) {
      tw_stime now = tw_now(lpptr);
      tw_stime trecv;

      double dt = tw_rand_exponential(lpptr->rng,mean_delay);
      trecv.t = now.t + dt;
      trecv.bits[0] = lpptr->gid;
      assert(trecv > now);
      if(trecv.t < g_tw_ts_end) {
	tw_event *evtptr = tw_event_new(lpptr->gid,trecv,lpptr);
	Event *evt = static_cast<Event *>(tw_event_data(evtptr));
	matcopy(n,evt->A,T);
	if(debug)
	  assert(invcheck(n,evt->A));
	tw_event_send(evtptr);
      }
    }
  }

  ~matlp() {
    global_stat.nforward += nforward;
    global_stat.nbackward += nbackward;
    global_stat.ncommit += ncommit;


    {
      uint tmp[n*n];
      matmul(n,global_stat.A,M,tmp);
      matcopy(n,global_stat.A,tmp);
    }

    delete[] T;
    delete[] M;
  }

  // typedef void (*init_f) (void *sv, tw_lp * me);
  static void init(void *state,tw_lp *twlp) {
    (void) new(state) matlp(twlp,nmat);
  }
  static void finalize(void *state,tw_lp *twlp) {
    static_cast<lp *>(state)->~lp();
  }

  void forward(tw_lp *twlp,Event *evt,const Time& now) {
    apply(evt->A);
    nforward++;

    /* Make new event to send... */ {
      long long int nlp_total =
        ((long long int) tw_nnodes()) *
        ((long long int) g_tw_npe) *
        ((long long int) g_tw_nlp);
      tw_stime now = tw_now(twlp);
      tw_stime trecv;
      
      double dt = tw_rand_exponential(twlp->rng,mean_delay);
      tw_lpid dest = tw_rand_integer(twlp->rng,0,nlp_total-1);
      trecv.t = now.t + dt;
      trecv.bits[0] = twlp->gid;
      assert(trecv > now);
      if(trecv.t < g_tw_ts_end) {
        tw_event *evtptr = tw_event_new(dest,trecv,twlp);
        Event *evt = static_cast<Event *>(tw_event_data(evtptr));
	matcopy(n,evt->A,T);
	if(debug)
	  assert(invcheck(n,evt->A));
        tw_event_send(evtptr);
      }
    }
  }
  void backward(tw_lp *twlp,Event *evt,const Time& now) {
    tw_rand_reverse_unif(twlp->rng);
    tw_rand_reverse_unif(twlp->rng);
    
    unapply(evt->A);
    nbackward++;
  }
  void commit(tw_lp *twlp,Event *evt,const Time& now) {
    ncommit++;
  }
};

#endif
