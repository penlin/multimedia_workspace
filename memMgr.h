#ifndef MEM_MGR_H
#define MEM_MGR_H

#include <cassert>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "build_value.h"

using namespace std;

// Turn this on for debugging
// #define MEM_DEBUG

#define SIZE_T      sizeof(size_t)
#define SIZE_T_1    (sizeof(size_t) - 1)

#define toSizeT(t)      (((t)%SIZE_T)?(t)-(t)%SIZE_T+SIZE_T:(t)-(t)%SIZE_T)
//
// To demote 't' to the nearest multiple of SIZE_T
// e.g. Let SIZE_T = 8;  downtoSizeT(9) = 8, downtoSizeT(100) = 96
#define downtoSizeT(t)  (t)-(t)%SIZE_T

// R_SIZE is the size of the recycle list
#define R_SIZE 64

//--------------------------------------------------------------------------
// Forward declarations
//--------------------------------------------------------------------------
class MemMgr;


//--------------------------------------------------------------------------
// Class Definitions
//--------------------------------------------------------------------------
// T is the class that use this memory manager
//
// Make it a private class;
// Only friend to MemMgr;
//
class MemBlock
{
   friend class MemMgr;

   // Constructor/Destructor
   MemBlock(MemBlock* n, size_t b) : _nextBlock(n) {
      _begin = _ptr = (char*)malloc(sizeof(char)*b);
      _end = _begin + b;
      memset(_begin,0,b);
  }
   ~MemBlock() {
        free(_begin);
   }

   // Member functions
   void reset() { _ptr = _begin; }

   bool getMem(size_t t, void*& ret) {
        size_t nt = toSizeT(t);
        if(nt<=getRemainSize()){
            ret= (void*) _ptr;
            _ptr= _ptr+nt;
            return true;
        } else {
            ret= (void*) _ptr;
            return false;
        }

   }
   size_t getRemainSize() const { return size_t(_end - _ptr); }

   MemBlock* getNextBlock() const { return _nextBlock; }

   // Data members
   char*             _begin;
   char*             _ptr;
   char*             _end;
   MemBlock*         _nextBlock;
};

// Make it a private class;
// Only friend to MemMgr;
//
class MemRecycleList
{
   friend class MemMgr;

   // Constructor/Destructor
   MemRecycleList(size_t a = 0) : _arrSize(a), _first(0), _nextList(0) {}
   ~MemRecycleList() { reset(); }

   // Member functions
   // ----------------
   // pop out the first element in the recycle list
   void* popFront() {
      // TODO
        if(_first!=0){
            void* temp = _first;
            _first=getNext(_first);
            return temp;
        }
        return 0;
   }
   // push the element 'p' to the beginning of the recycle list
    void  pushFront(void* p) {
      // TODO
        if(_first!=0){
           *((size_t *)p) = ((size_t)_first);
            _first = p;
        }
        else {
            _first = p;
            *((size_t *)p)=0;
        }
    }
   // Release the memory occupied by the recycle list(s)
   // DO NOT release the memory occupied by MemMgr/MemBlock
   void reset() {
        _first =0;
        if(_nextList != 0) _nextList->reset();
   }

   // Helper functions
   // ----------------
   // Get the next element after 'p' in the recycle list
   void* getNext(void* p) const {
      // TODO
        if((*(size_t *)p)!=0)
            return (void *)(*(size_t *)p);
        else return 0;
   }

   MemRecycleList* getList(size_t n) {
        assert((n - _arrSize%R_SIZE) % R_SIZE == 0);
        if(_arrSize==n) return this;
        else if(_nextList!=0) return _nextList->getList(n);
        else if(_first==0){
            _arrSize = n;
            return this;
        }
        else {
            _nextList=new MemRecycleList(n);
            return _nextList;
        }
   }
   // count the number of elements in the recycle list
   size_t numElm() const {
        size_t count = 0;
        void* p = _first;
        while (p) {
             p = getNext(p);
             ++count;
        }
        return count;
   }

   // Data members
   size_t              _arrSize;   // the array size of the recycled data
   void*               _first;     // the first recycled data
   MemRecycleList*     _nextList;  // next MemRecycleList
                                   //      with _arrSize + x*R_SIZE
};

class MemMgr
{

public:
   MemMgr(size_t b = 65536) : _blockSize(b) {
      assert(b % SIZE_T == 0);
      _activeBlock = new MemBlock(0, _blockSize);
	  int i=0;
      for (i = 0; i < R_SIZE; ++i)
         _recycleList[i]._arrSize = i;
#if !__MEM_MGR_BOOST__
      _max_alloc = 0;
      _min_alloc = b;
#endif
   }
   ~MemMgr() { reset(); delete _activeBlock; }

   // 1. Remove the memory of all but the firstly allocated MemBlocks
   //    That is, the last MemBlock searchd from _activeBlock.
   //    reset its _ptr = _begin (by calling MemBlock::reset())
   // 2. reset _recycleList[]
   // 3. 'b' is the new _blockSize; "b = 0" means _blockSize does not change
   //    if (b != _blockSize) reallocate the memory for the first MemBlock
   // 4. Update the _activeBlock pointer
   void reset(size_t b = 0) {
#if MEM_DEBUG
        cout << "Resetting memMgr...(" << b << ")" << endl;
#endif // MEM_DEBUG
        while(_activeBlock->getNextBlock()!=NULL){
            MemBlock* temp =_activeBlock;
            _activeBlock=_activeBlock->getNextBlock();
            delete temp;
        }
        _activeBlock->reset();
        int i=0;
        for(i=0;i<R_SIZE;i++){
            _recycleList[i].reset();
        }
        if(b!=_blockSize && b!=0) {
            _blockSize=b;
            _activeBlock = new MemBlock(0, _blockSize);
        }
   }
   // Called by new
   void* alloc(const size_t &t) {
//        #ifdef MEM_DEBUG
//        cout << "Calling alloc...(" << t+SIZE_T << ")" << endl;
//        #endif // MEM_DEBUG
#if !__MEM_MGR_BOOST__
        if(!t%SIZE_T)
            printf("want to alloc size %d\n",t);
#endif
        void* ret = getMem(t+SIZE_T);
        *(size_t *)ret = (t+SIZE_T);
        return ((size_t *)ret )+1;
   }
   // Called by new[]
   void* allocArr(const size_t &t, const size_t &n) {
//        #ifdef MEM_DEBUG
//        cout << "Calling allocArr...(" << t+SIZE_T << ")" << endl;
//        #endif // MEM_DEBUG
        // Note: no need to record the size of the array == > system will do
        return alloc(n*t);
   }

   void* alloc2DMat(const size_t & t, const size_t & r, const size_t & c ){
        size_t* ret = (size_t *)alloc(r*c*t+r*SIZE_T);
//        ret[0] = (size_t)(ret+r);
//        for(int i = 0, step = c*t ; i < r ; ++i)
//            printf("ret[%d]=%p  , &ret[%d]=%p\n",i,ret[i],i+r,&ret[i+r]);
        return (void*) ret;
   }

   void* alloc3DMat(const size_t & t, const size_t & r, const size_t & c, const size_t &d ){
        size_t sec_len = r*(c+1);
        size_t* ret = (size_t *)alloc(r*c*d*t+sec_len*SIZE_T);
//        int i = 1, step = c*SIZE_T;
//        ret[0] = (size_t)(ret+r);
//        for(; i < r ; ++i)
//            ret[i] = (size_t)(((char*)ret[i-1])+step);
//        ret[r] = (size_t)(ret+sec_len);
//        for( i = r+1, step = d*t ; i < sec_len ; ++ i)
//            ret[i] = (size_t)(((char*)ret[i-1])+step);
        return (void*) ret;
   }
   // Called by delete
   void  free(void* p) {
        void* ret = (void *)(((size_t *)p)-1);
        size_t n =*(size_t*)ret;
        #if MEM_DEBUG
        cout << "Calling free...(" << n << ")" << endl;
        #endif // MEM_DEBUG
        getMemRecycleList(getRecycleIdx(n))->pushFront(ret);
   }
//   // Called by delete[]
//   void  freeArr(void* p) {
//        void* ret = (void *)(((size_t *)p)-1);
//        size_t n =*(size_t*)ret;
//        #ifdef MEM_DEBUG
//        cout << "Calling free...(" << n << ")" << endl;
//        #endif // MEM_DEBUG
//        getMemRecycleList(getRecycleIdx(n))->pushFront(ret);
//      // TODO
//      // Get the array size 'n' stored by system,
//      // which is also the _recycleList index
//      // ==> assert(n == getRecycleIdx(n * S + SIZE_T));
//   }
   void print() const {
      cout << "=========================================" << endl
           << "=              Memory Manager           =" << endl
           << "=========================================" << endl
           << "* Block size            : " << _blockSize << " Bytes" << endl
           << "* Number of blocks      : " << getNumBlocks() << endl
           << "* Free mem in last block: " << _activeBlock->getRemainSize() << endl;
      int total = _blockSize*getNumBlocks() - _activeBlock->getRemainSize() ;
      cout << "* Total mem alloced     : " << total << " Bytes" << endl
#if !__MEM_MGR_BOOST__
           << "* Maximun/Minimum Memory Allocation: " << _max_alloc << " / " << _min_alloc << endl
#endif
           << endl
           << "* Recycle list          : " << endl;
      int i = 0, count = 0, recycle = 0;
      while (i < R_SIZE) {
         const MemRecycleList* ll = &(_recycleList[i]);
         while (ll != 0) {
            size_t s = ll->numElm();
            if (s) {
               cout << "[" << setw(3) << right << i << "(" << ll->_arrSize << ")] = "
                    << setw(10) << left << s;
               if (++count % 4 == 0) cout << endl;
               recycle+=(ll->_arrSize+1)*SIZE_T;
            }
            ll = ll->_nextList;
         }
         ++i;
      }
      cout << endl
           << "* Recycled Memory       : " << recycle << " Bytes" << endl
           << "* Inused Memory         : " << total-recycle << " Bytes" << endl;
   }

private:
   size_t                     _blockSize;
   MemBlock*                  _activeBlock;
   MemRecycleList             _recycleList[R_SIZE];
#if !__MEM_MGR_BOOST__
   size_t                     _max_alloc;
   size_t                     _min_alloc;
#endif
   // Private member functions
   //
   // t: #Bytes; MUST be a multiple of SIZE_T
   // return index for recycle list
   // [Note] t must >= S
   // [NOTE] Use this function in (at least) getMem() for retrieving/storing
   //        memory from recycle list
   size_t getRecycleIdx(size_t t) const {
        assert(t % SIZE_T == 0);
        return (t/SIZE_T-1);
   }
   // Once the index for recycle list is known, use this function to get the
   // recycle list whose _arrSize == n
   MemRecycleList* getMemRecycleList(size_t n) {
        size_t m = n % R_SIZE;
        return _recycleList[m].getList(n);
   }
   // t is the #Bytes requested from new or new[]
   // Note: Make sure the returned memory is a multiple of SIZE_T
   void* getMem(size_t t) {
        void* ret = 0;
        #if MEM_DEBUG
        cout << "Calling MemMgr::getMem...(" << t << ")" << endl;
        #endif // MEM_DEBUG
        size_t nt = toSizeT(t);
#if !__MEM_MGR_BOOST__
        if(nt>_max_alloc)
            _max_alloc = nt;
        if (nt < _min_alloc)
            _min_alloc = nt;
#endif
        if(nt>_blockSize){
            cerr << "Requested memory (" << t << ") is greater than block size"
             << "(" << _blockSize << "). " << "Exception raised...\n";
            throw bad_alloc();
        }

        MemRecycleList* temp = getMemRecycleList(getRecycleIdx(nt));
        if(temp->_first != 0){
            #if MEM_DEBUG
            printf("getMem from RecycleList[%d]\n",getRecycleIdx(nt));
            #endif // MEM_DEBUG
            ret = temp->popFront();
            return ret;
        }else if(_activeBlock->getMem(nt,ret)){
             #if MEM_DEBUG
             printf("getMem from Memmory Block\n");
             #endif // MEM_DEBUG
             return ret;
        }else {
        ///// *** recycle the remained memory *** ///////
            size_t remain = _activeBlock->getRemainSize();
            size_t recycle = downtoSizeT(remain);
            if(remain < 1+SIZE_T && remain >= 1) getMemRecycleList(0)->pushFront(ret);
            else if(remain >= 1+SIZE_T)
                //*(ret+recycle)
                getMemRecycleList(getRecycleIdx(recycle))->pushFront(ret);

        /////**** allocate a new memory block *****//////
            _activeBlock= new MemBlock(_activeBlock,_blockSize);
            if(_activeBlock->getMem(nt,ret))
                return ret;
        }

        return ret;
   }
   // Get the currently allocated number of MemBlock's
   size_t getNumBlocks() const {
      // TODO
        MemBlock* temp=_activeBlock;
        size_t count=1;
        while(temp->getNextBlock()!=NULL){
            count++;
            temp=temp->getNextBlock();
        }
        return count;
   }

};

#endif // MEM_MGR_H
