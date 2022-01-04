/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package btllib;

public class VectorString extends java.util.AbstractList<String> implements java.util.RandomAccess {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected VectorString(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(VectorString obj) {
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
        btllibJNI.delete_VectorString(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public VectorString(String[] initialElements) {
    this();
    reserve(initialElements.length);

    for (String element : initialElements) {
      add(element);
    }
  }

  public VectorString(Iterable<String> initialElements) {
    this();
    for (String element : initialElements) {
      add(element);
    }
  }

  public String get(int index) {
    return doGet(index);
  }

  public String set(int index, String e) {
    return doSet(index, e);
  }

  public boolean add(String e) {
    modCount++;
    doAdd(e);
    return true;
  }

  public void add(int index, String e) {
    modCount++;
    doAdd(index, e);
  }

  public String remove(int index) {
    modCount++;
    return doRemove(index);
  }

  protected void removeRange(int fromIndex, int toIndex) {
    modCount++;
    doRemoveRange(fromIndex, toIndex);
  }

  public int size() {
    return doSize();
  }

  public VectorString() {
    this(btllibJNI.new_VectorString__SWIG_0(), true);
  }

  public VectorString(VectorString other) {
    this(btllibJNI.new_VectorString__SWIG_1(VectorString.getCPtr(other), other), true);
  }

  public long capacity() {
    return btllibJNI.VectorString_capacity(swigCPtr, this);
  }

  public void reserve(long n) {
    btllibJNI.VectorString_reserve(swigCPtr, this, n);
  }

  public boolean isEmpty() {
    return btllibJNI.VectorString_isEmpty(swigCPtr, this);
  }

  public void clear() {
    btllibJNI.VectorString_clear(swigCPtr, this);
  }

  public VectorString(int count, String value) {
    this(btllibJNI.new_VectorString__SWIG_2(count, value), true);
  }

  private int doSize() {
    return btllibJNI.VectorString_doSize(swigCPtr, this);
  }

  private void doAdd(String x) {
    btllibJNI.VectorString_doAdd__SWIG_0(swigCPtr, this, x);
  }

  private void doAdd(int index, String x) {
    btllibJNI.VectorString_doAdd__SWIG_1(swigCPtr, this, index, x);
  }

  private String doRemove(int index) {
    return btllibJNI.VectorString_doRemove(swigCPtr, this, index);
  }

  private String doGet(int index) {
    return btllibJNI.VectorString_doGet(swigCPtr, this, index);
  }

  private String doSet(int index, String val) {
    return btllibJNI.VectorString_doSet(swigCPtr, this, index, val);
  }

  private void doRemoveRange(int fromIndex, int toIndex) {
    btllibJNI.VectorString_doRemoveRange(swigCPtr, this, fromIndex, toIndex);
  }

}
