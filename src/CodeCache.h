/*
 * Copyright (c) Facebook, Inc. and its affiliates.
 * All rights reserved.
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#pragma once
#include <map>

namespace fbgemm {

/**
 * @brief Thread safe cache for microkernels, ensures single creation per key.
 * @tparam Key Type of unique key (typically a tuple)
 * @tparam Value Type of the microkernel function (Typically a function pointer)
 */
template <typename KEY, typename VALUE>
class CodeCache {
 private:
  std::map<KEY, VALUE> values_;

 public:
  CodeCache(const CodeCache&) = delete;
  CodeCache& operator=(const CodeCache&) = delete;

  CodeCache(){};

  VALUE getOrCreate(const KEY& key, std::function<VALUE()> generatorFunction) {
    auto it = values_.find(key);
    if (it != values_.end()) {
      return it->second;
    } else {
      //std::cerr << "create" << std::endl;
      auto fn = generatorFunction();
      values_[key] = fn;
      return fn;
    }
  }
};

} // namespace fbgemm