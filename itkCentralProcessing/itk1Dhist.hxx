#ifndef ITK1DHIST_HEADER
#define ITK1DHIST_HEADER

#include <vector>

template <typename Tdatain, typename pdfT, int Nbins> class hist1D {
protected:
public:
  typedef unsigned int HISTT;
  std::vector<HISTT> hist;
  std::vector<pdfT> prob;
  Tdatain min,max,range,binsize;
  typedef std::size_t BINIDX;
  hist1D(void) : hist(Nbins),prob(Nbins),min(0),max(0),range(0),binsize(0) {
    clear();
  }

  void initialize(Tdatain _min, Tdatain _max) {
    min=_min;
    max=_max;
    range=max-min;
    binsize=range/Nbins;
  }

  void clear(void) {
    for (int i=0;i<Nbins;i++) {
      hist[i]=0;
      prob[i]=0;
    }
  }

  void add(Tdatain value) {
    assert(max>min);
    BINIDX b=bin(value);
    hist[b]++;
  }

  void normalize(void) {
    assert(max>min);
    pdfT total=0;
    for (int i=0;i<Nbins;i++) 
      total+=hist[i];
    if (total>0) {
      for (int i=0;i<Nbins;i++) 
	prob[i]=hist[i]/total;
    }
  }

  BINIDX bin(Tdatain value) {
    BINIDX b=0;
    if (value>=max)
      b=Nbins-1;
    else if (value<=min)
      b=0;
    else 
      b=BINIDX((value-min)/binsize);
    return b;
  }

  pdfT operator [] (BINIDX i) const { return prob[i]; } //get

  Tdatain midbinvalue(BINIDX i) const { return i*binsize+min+binsize/2.; }
};


template <typename Tdatain, int Qdims, typename pdfT, int Nbins> class histQx1D {
protected:
public:
  typedef hist1D<Tdatain, pdfT, Nbins> HISTST;
  std::vector<HISTST> hists;
  typedef std::size_t BINIDX;
  histQx1D(void): hists(Qdims) {
  }
  template <typename structInT> void add(const structInT &values) {
    for (int i=0;i<Qdims;i++)
      hists[i].add(values[i]);
  }
  void normalize(void) {
    for (int i=0;i<Qdims;i++)
      hists[i].normalize();
  }
  void clear(void) {
    for (int i=0;i<Qdims;i++)
      hists[i].clear();
  }
  void initialize(Tdatain _min, Tdatain _max) {
    for (int i=0;i<Qdims;i++)
      hists[i].initialize(_min,_max);
  }

  HISTST operator [] (int i) const { return hists[i]; }
  HISTST & operator [] (int i) { return hists[i]; }
};

#endif //ITK1DHIST_HEADER
