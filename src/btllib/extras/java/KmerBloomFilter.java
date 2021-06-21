/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package btllib;

public class KmerBloomFilter {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected KmerBloomFilter(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(KmerBloomFilter obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  @SuppressWarnings("deprecation")
  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        btllibJNI.delete_KmerBloomFilter(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public KmerBloomFilter() {
    this(btllibJNI.new_KmerBloomFilter__SWIG_0(), true);
  }

  public KmerBloomFilter(long bytes, long hash_num, long k) {
    this(btllibJNI.new_KmerBloomFilter__SWIG_1(bytes, hash_num, k), true);
  }

  public KmerBloomFilter(String path) {
    this(btllibJNI.new_KmerBloomFilter__SWIG_2(path), true);
  }

  public void insert(String seq, long seq_len) {
    btllibJNI.KmerBloomFilter_insert__SWIG_0(swigCPtr, this, seq, seq_len);
  }

  public void insert(String seq) {
    btllibJNI.KmerBloomFilter_insert__SWIG_1(swigCPtr, this, seq);
  }

  public void insert(SWIGTYPE_p_uint64_t hashes) {
    btllibJNI.KmerBloomFilter_insert__SWIG_2(swigCPtr, this, SWIGTYPE_p_uint64_t.getCPtr(hashes));
  }

  public void insert(SWIGTYPE_p_std__vectorT_uint64_t_t hashes) {
    btllibJNI.KmerBloomFilter_insert__SWIG_3(swigCPtr, this, SWIGTYPE_p_std__vectorT_uint64_t_t.getCPtr(hashes));
  }

  public long contains(String seq, long seq_len) {
    return btllibJNI.KmerBloomFilter_contains__SWIG_0(swigCPtr, this, seq, seq_len);
  }

  public long contains(String seq) {
    return btllibJNI.KmerBloomFilter_contains__SWIG_1(swigCPtr, this, seq);
  }

  public boolean contains(SWIGTYPE_p_uint64_t hashes) {
    return btllibJNI.KmerBloomFilter_contains__SWIG_2(swigCPtr, this, SWIGTYPE_p_uint64_t.getCPtr(hashes));
  }

  public boolean contains(SWIGTYPE_p_std__vectorT_uint64_t_t hashes) {
    return btllibJNI.KmerBloomFilter_contains__SWIG_3(swigCPtr, this, SWIGTYPE_p_std__vectorT_uint64_t_t.getCPtr(hashes));
  }

  public long get_bytes() {
    return btllibJNI.KmerBloomFilter_get_bytes(swigCPtr, this);
  }

  public SWIGTYPE_p_uint64_t get_pop_cnt() {
    return new SWIGTYPE_p_uint64_t(btllibJNI.KmerBloomFilter_get_pop_cnt(swigCPtr, this), true);
  }

  public double get_occupancy() {
    return btllibJNI.KmerBloomFilter_get_occupancy(swigCPtr, this);
  }

  public long get_hash_num() {
    return btllibJNI.KmerBloomFilter_get_hash_num(swigCPtr, this);
  }

  public double get_fpr() {
    return btllibJNI.KmerBloomFilter_get_fpr(swigCPtr, this);
  }

  public long get_k() {
    return btllibJNI.KmerBloomFilter_get_k(swigCPtr, this);
  }

  public String get_hash_fn() {
    return btllibJNI.KmerBloomFilter_get_hash_fn(swigCPtr, this);
  }

  public BloomFilter get_bloom_filter() {
    return new BloomFilter(btllibJNI.KmerBloomFilter_get_bloom_filter(swigCPtr, this), false);
  }

  public void save(String path) {
    btllibJNI.KmerBloomFilter_save(swigCPtr, this, path);
  }

}
