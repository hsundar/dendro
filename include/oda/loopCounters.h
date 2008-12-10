
#ifndef __LOOP_COUNTERS_H__
#define __LOOP_COUNTERS_H__

#include "Point.h"

namespace ot {
        struct LoopCounters {
          Point currentOffset;
          unsigned int currentIndex;
          unsigned int qCounter;
          unsigned int pgQcounter;
        };
} //end namespace

#endif

