/*
 * Copyright (c) Facebook, Inc. and its affiliates.
 * All rights reserved.
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#pragma once
#include <condition_variable>
#include <future>
#include <unordered_map>

#if __cplusplus >= 201402L && !defined(__APPLE__)
// For C++14, use shared_timed_mutex.
// some macOS C++14 compilers don't support shared_timed_mutex.
#define FBGEMM_USE_SHARED_TIMED_MUTEX
#endif

#ifdef FBGEMM_USE_SHARED_TIMED_MUTEX
#include <shared_mutex>
#else
#include <mutex>
#endif

namespace fbgemm {

template <class T> using hash = std::hash<T>;

// This combinator is based on boost::hash_combine, but uses
// std::hash as the hash implementation. Used as a drop-in
// replacement for boost::hash_combine.

template <class T>
inline void hash_combine(std::size_t& seed, T const& v) {
  hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/**
 * @brief Thread safe cache for microkernels, ensures single creation per key.
 * @tparam Key Type of unique key (typically a tuple)
 * @tparam Value Type of the microkernel function (Typically a function pointer)
 */
template <typename KEY, typename VALUE>
class CodeCache {
 private:
  std::unordered_map<KEY, std::shared_future<VALUE>> values_;
#ifdef FBGEMM_USE_SHARED_TIMED_MUTEX
  std::shared_timed_mutex mutex_;
#else
  std::mutex mutex_;
#endif

 public:
  CodeCache(const CodeCache&) = delete;
  CodeCache& operator=(const CodeCache&) = delete;

  CodeCache(){};

  VALUE getOrCreate(const KEY& key, std::function<VALUE()> generatorFunction) {
    std::shared_future<VALUE> returnFuture;
    std::promise<VALUE> returnPromise;
    bool needsToGenerate = false;

    // Check for existance of the key
    {
#ifdef FBGEMM_USE_SHARED_TIMED_MUTEX
      mutex_.lock_shared();
#else
      std::unique_lock<std::mutex> lock(mutex_);
#endif

      auto it = values_.find(key);
      if (it != values_.end()) {
        returnFuture = it->second;
#ifdef FBGEMM_USE_SHARED_TIMED_MUTEX
        mutex_.unlock_shared();
#endif
      } else {
#ifdef FBGEMM_USE_SHARED_TIMED_MUTEX
        mutex_.unlock_shared();

        mutex_.lock();
        // Need to look up again because there could be race condition from
        // the time gap between mutex_.unlock_shared() and mutex_.lock()
        it = values_.find(key);
        if (it == values_.end()) {
#endif
          values_[key] = returnFuture = returnPromise.get_future().share();
          needsToGenerate = true;
#ifdef FBGEMM_USE_SHARED_TIMED_MUTEX
        } else {
          returnFuture = it->second;
        }
        mutex_.unlock();
#endif
      }
    }

    // The value (code) generation is not happening under a lock
    if (needsToGenerate) {
      returnPromise.set_value(generatorFunction());
    }

    // Wait for the future and return the value
    return returnFuture.get();
  }
};

} // namespace fbgemm