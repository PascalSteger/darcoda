//MySTL (c) 2008 by Oliver Hahn, hahn@phys.ethz.ch

/*!
 * \file MySTL.h
 * \brief Customized STL functions
 */

#ifndef MYSTL_H
#define MYSTL_H

#include <cstdlib>
#include <list>
#include <iostream>

typedef unsigned Identifier;

template<typename T>
class array{
  friend class GadgetSnapshot;
 protected:
  T *m_Data;
  Identifier m_Size;
  Identifier m_Backinsert;
 public:
  array( void ) : m_Data(NULL), m_Size(0), m_Backinsert(0) {}; 
  array( Identifier aSize ){
    m_Data = new T[aSize];
    m_Size = aSize;
    m_Backinsert = 0;
  }
  array( array<T> &aArray ){
    m_Size = aArray.m_Size;
    m_Data = new T[m_Size];
    m_Backinsert = aArray.m_Size;
    for( Identifier i=0; i<m_Size; ++i )
      m_Data[i] = aArray[i];
  }
  ~array(){
    if(m_Data != NULL)
      delete[] m_Data;
  }
  // just encapsulate a given array
  array( T* dataptr, Identifier size ){
    m_Data = dataptr;
    m_Size = size;
  }
  inline void reserve( Identifier aSize ){
    if( m_Data != NULL )
      delete[] m_Data;
    m_Size = aSize;
    m_Data = new T[aSize];
  }
  inline void resize( Identifier aSize ){
    T *DataTemp = new T[aSize];
    for( Identifier i=0; i<m_Size && i<aSize; ++i )
      DataTemp[i] = m_Data[i];
    delete[] m_Data;
    m_Data = DataTemp;
    m_Size = aSize;
  }
		
  inline void push_back( T Elem ){
    m_Data[m_Backinsert++] = Elem;
  }
  inline T &operator[]( Identifier i ) const{
    return m_Data[i];
  }
  inline Identifier size( void ) const{
    return m_Size;
  }
  inline bool find( T elem ) const{
    bool bRetVal = false;
    for( Identifier i=0; i<m_Size; ++i ){
      if( elem == m_Data[i] ){
	bRetVal = true;
	break;
      }
    }
    return bRetVal;
  }
  inline void clear( void ) {
    m_Size = 0;
    delete[] m_Data;
    m_Data = NULL;
  }
  inline T* data_ptr( void ){
    return &m_Data[0];
  }
};



#include <iterator>
template <class T>
class list
{
  struct Node
  {
    Node(const T& x,Node* y = 0):m_data(x),m_next(y){}
    T m_data;
    Node* m_next;
  };

  Node* m_node;
  Identifier m_n;

 public:
  // need forward declaration to avoid warnings
  class const_iterator;

  class iterator
    : public std::iterator<std::forward_iterator_tag, T>
    {
      Node* m_rep;
    public:
      friend class list;
      friend class const_iterator;
      typedef T value_type;
      typedef T& reference;
      typedef T* pointer;
      typedef int difference_type;
      typedef std::forward_iterator_tag iterator_category;

      inline iterator(Node* x=0):m_rep(x){}
      inline iterator(const iterator& x):m_rep(x.m_rep) {}
      inline iterator& operator=(const iterator& x)
	{ 
	  m_rep=x.m_rep; return *this; 
	}
      inline iterator& operator++()
	{ 
	  m_rep = m_rep->m_next; return *this; 
	}
      inline iterator operator++(int)
	{ 
	  iterator tmp(*this); m_rep = m_rep->m_next; return tmp; 
	}
      inline reference operator*() const { return m_rep->m_next->m_data; }
      inline pointer operator->() const { return m_rep->m_next; }
      inline bool operator==(const iterator& x) const
	{
	  return m_rep == x.m_rep; 
	}	
      inline bool operator!=(const iterator& x) const
	{
	  return m_rep != x.m_rep; 
	}	

    };

  class const_iterator
    : public std::iterator<std::forward_iterator_tag, const T>
    {
      const Node* m_rep;
    public:
      friend class list;
      friend class iterator;
      typedef T value_type;
      typedef T& reference;
      typedef T* pointer;
      inline const_iterator(const Node* x=0):m_rep(x){}
      inline const_iterator(const const_iterator& x):m_rep(x.m_rep) {}
      inline const_iterator(const iterator& x):m_rep(x.m_rep){}
      inline const_iterator& operator=(const const_iterator& x)
	{ 
	  m_rep=x.m_rep; return *this; 
	}
      inline const_iterator& operator=(const iterator& x)
	{ 
	  m_rep=x.m_rep; return *this; 
	}		
      inline const_iterator& operator++()
	{ 
	  m_rep = m_rep->m_next; return *this; 
	}
      inline const_iterator operator++(int)
	{ 
	  const_iterator tmp(*this); m_rep = m_rep->m_next; return tmp; 
	}
      inline reference operator*() const { return m_rep->m_next->m_data; }
      inline pointer operator->() const { return m_rep->m_next; }
      inline bool operator==(const const_iterator& x) const
	{
	  return m_rep == x.m_rep; 
	}
      inline bool operator!=(const const_iterator& x) const
	{
	  return m_rep != x.m_rep; 
	}



    };


  list() : m_node(new Node(T())) { m_node->m_next = m_node; m_n=0;}

  list(const list& L) : m_node(new Node(T()))
    {
      m_n = 0;
      m_node->m_next=m_node;
      for ( const_iterator i = L.begin(); i!= L.end(); ++i )
	push_front(*i);
      reverse();
    }

  void reverse()
    {
      if (empty())
	return;
      Node* new_m_node = m_node->m_next->m_next;
      Node* p = m_node; Node* i = m_node->m_next; Node* n;
      do  
	{
	  n = i->m_next;
	  i->m_next = p;
	  p = i; i = n;
	}
      while (p!=m_node);
      m_node = new_m_node;
    }

  void swap(list& x)
    {
      Node* tmp = m_node; m_node = x.m_node; x.m_node = tmp;
      Identifier ntmp = m_n; m_n = x.m_n; x.m_n = ntmp;
    }

  list& operator=(const list& x)
    {
      list tmp(x);
      swap(tmp);
      return *this;
    }

  ~list() { clear(); delete m_node; }
  //	void clear() { while (!empty()) pop_front(); }
  inline void clear() {
    Node* cur = (Node*) m_node->m_next;
    while( cur != m_node ){
      Node *tmp = cur;
      cur = (Node*) cur->m_next;
      //destroy( &tmp->m_data );
      //put_node( tmp );
      delete tmp;
    }
    m_node->m_next = m_node;
    m_n = 0;
  }
			



  inline void push_front(const T&x)
    {
      insert (begin(),x);
    }
  inline void push_back(const T& x)
    {
      insert (end(),x);
    }
  inline void pop_front()
    {
      erase(begin());
    }
  //inline bool empty() { return m_node == m_node->m_next; }
  inline bool empty() { return m_n == 0; }

  inline T& front() { return *begin(); }
  inline const T& front() const { return *begin(); }
  inline T& back() { return m_node->data; }
  inline const T& back() const { return m_node->data; }

  inline iterator begin() { return iterator(m_node->m_next); }
  //inline iterator begin() { return (Node*)(m_node->m_next); }
  inline iterator end() { return iterator(m_node); }
  //inline iterator end() { return (Node*)m_node; }
  inline const_iterator begin() const 
    { 
      return const_iterator(m_node->m_next); 
    }
  inline const_iterator end() const 
    { 
      return const_iterator(m_node); 
    }

  iterator erase (iterator x)
    {
      if (x==end())
	return x;
      if (x.m_rep->m_next == m_node)
	m_node = x.m_rep;

      Node* tmp = x.m_rep->m_next;
      x.m_rep->m_next = x.m_rep->m_next->m_next;
      delete tmp;
      --m_n;
      return x;
    }

  iterator insert (iterator x, const T& y)
    {
      Node* tmp = new Node(y,x.m_rep->m_next);
      if (x.m_rep == m_node)
	m_node = tmp;
      x.m_rep->m_next = tmp;
      ++m_n;
      return x;
    }

  // rotate x to beginning
  void rotate (iterator x)
    {
      if (x == end())
	return;
      Node* sentinel = m_node->m_next;
      m_node->m_next = m_node->m_next->m_next;
      sentinel->m_next = x.m_rep->m_next;
      x.m_rep->m_next = sentinel;
      m_node = x.m_rep; 
    }
  /*void splice( list<T>& other )
    {
    iterator x=other.begin();
    while( x!=other.end() ){
    push_back( *x );
    other.erase(x);
    }
		
    //end().m_rep->m_next = other.m_node->m_next;
    //other.m_node->m_next = other.m_node;
		
		
    //		other.clear();
    }*/
	
  inline void splice( list<T>& other )
    {		
      Node *head1, *head2;
      head1 = m_node->m_next;
      head2 = other.m_node->m_next;
		
      m_node->m_next = other.m_node->m_next->m_next;
      other.m_node->m_next = head1;
      m_node = other.m_node;
		
      other.m_node = head2;
      head2->m_next = head2;
		
      m_n += other.m_n;
      other.m_n = 0;
    }
	
  inline Identifier size( void )
    {
      return m_n;
    }
			
};

template<class T>
void 
convertlists( std::list<T> &listto, list<T> &listfrom )
{
  typename list<T>::iterator it = listfrom.begin();
  while( it!=listfrom.end() ){
    listto.push_back(*it);
    ++it;
  }
}

#endif /*MYSTL_HH_*/
