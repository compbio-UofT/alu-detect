#ifndef concurrent_queue_hpp_
#define concurrent_queue_hpp_

#include <iostream>
#include <queue>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>


template<typename Data>
class concurrent_queue
{
private:
  static int const initial_capacity = (1 << 10);
  static int const capacity_increase = (1 << 10);

  std::queue<Data> the_queue;
  mutable boost::mutex the_mutex;
  boost::condition_variable the_condition_variable;
  bool done;
  size_t capacity;

  class queue_not_empty_or_done
  {
  public:
    std::queue<Data>& queue;
    bool& done;
    queue_not_empty_or_done(std::queue<Data>& queue_, bool& done_)
      : queue(queue_), done(done_) {}
    bool operator()() const { return not queue.empty() or done; }
  };

  class queue_not_full
  {
  public:
    std::queue<Data>& queue;
    size_t& capacity;
    queue_not_full(std::queue<Data>& queue_, size_t& capacity_)
      : queue(queue_), capacity(capacity_) {
    }
    bool operator()() const {
      return queue.size() < capacity;
    }
  };


public:
  concurrent_queue() : done(false), capacity(initial_capacity) {}

  int timed_wait_and_push(Data const& data, int s = 2)
  {
    boost::mutex::scoped_lock lock(the_mutex);
    boost::posix_time::time_duration td = boost::posix_time::seconds(s);
    the_condition_variable.timed_wait(lock, td, queue_not_full(the_queue, capacity));

    int ret = 0;
    if (the_queue.size() >= capacity) {
      capacity += capacity_increase;
      ret = 1;
    }
    the_queue.push(data);

    lock.unlock();
    the_condition_variable.notify_one();
    return ret;
  }

  void wait_and_push(Data const& data)
  {
    boost::mutex::scoped_lock lock(the_mutex);
    the_condition_variable.wait(lock, queue_not_full(the_queue, capacity));

    the_queue.push(data);

    lock.unlock();
    the_condition_variable.notify_one();
  }

  bool empty() const
  {
    boost::mutex::scoped_lock lock(the_mutex);
    return the_queue.empty();
  }

  bool try_pop(Data& popped_value)
  {
    boost::mutex::scoped_lock lock(the_mutex);
    if (the_queue.empty()) {
      return false;
    }
        
    popped_value = the_queue.front();
    the_queue.pop();
    lock.unlock();
    the_condition_variable.notify_one();
    return true;
  }

  bool wait_and_pop(Data& popped_value)
  {
    boost::mutex::scoped_lock lock(the_mutex);
    the_condition_variable.wait(lock, queue_not_empty_or_done(the_queue, done));

    bool res = false;
    if (not the_queue.empty()) {
      popped_value = the_queue.front();
      the_queue.pop();
      res = true;
    }

    lock.unlock();
    the_condition_variable.notify_one();
    return res;
  }

  void set_done()
  {
    boost::mutex::scoped_lock lock(the_mutex);
    done = true;
    lock.unlock();
    the_condition_variable.notify_one();
  }
};


#endif
