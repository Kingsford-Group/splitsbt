#ifndef __BIT_POINTER_H__
#define __BIT_POINTER_H__

#include <cstdint>
#include <iterator>

// A pointer to bits. Can store boolean (bits). Increment (++) advance
// by 1 bit only. It is marked as a random access iterator, although
// the implementation is not complete (TBD). This is the const version
// (read only pointer).
class const_bit_pointer
  : public std::iterator<std::random_access_iterator_tag, bool>
{
protected:
  uint64_t* ptr_m;
  uint64_t  mask_m;

public:
  const_bit_pointer(uint64_t* ptr)
    : ptr_m(ptr)
    , mask_m(1)
  { }

  const_bit_pointer& operator++() {
    mask_m <<= 1;
    if(!mask_m) {
      ++ptr_m;
      mask_m = 1;
    }
    return *this;
  }
  const_bit_pointer operator++(int) {
    const_bit_pointer res(*this);
    ++*this;
    return res;
  }
  bool operator==(const const_bit_pointer& rhs) const {
    return ptr_m == rhs.ptr_m && mask_m == rhs.mask_m;
  }
  bool operator<(const const_bit_pointer& rhs) const {
    return ptr_m < rhs.ptr_m || (ptr_m == rhs.ptr_m && mask_m < rhs.mask_m);
  }
  bool operator*() const { return (*ptr_m & mask_m) != 0; }

    uint64_t* get_ptr() const{
        return ptr_m;
    }

    uint64_t get_mask() const{
        return mask_m;
    }

};

// This is the read/write version of the bit pointer.
class bit_pointer
  : public const_bit_pointer
{
  class setter : public const_bit_pointer {
  public:
    setter(const bit_pointer& rhs)
      : const_bit_pointer((const_bit_pointer)rhs)
    { }
    setter& operator=(bool x) {
      if(x)
        *ptr_m |= mask_m;
      else
        *ptr_m &= ~mask_m;
      return *this;
    }
    operator bool() const { return const_bit_pointer::operator*(); }
    bit_pointer operator&() { return bit_pointer(*this); }
  };
public:
  bit_pointer(uint64_t* ptr)
    : const_bit_pointer(ptr)
  { }
  bit_pointer(const setter& rhs)
    : const_bit_pointer((const_bit_pointer)rhs)
  { }

  bit_pointer& operator++() {
    const_bit_pointer::operator++();
    return *this;
  }
  bit_pointer operator++(int) {
    bit_pointer res(*this);
    ++*this;
    return res;
  }
  setter operator*() { return setter(*this); }

    uint64_t* get_ptr() const{
        return ptr_m;
    }

    uint64_t get_mask() const{
        return mask_m;
    }

};



#endif /* __BIT_POINTER_H__ */
