
/**
  @file dendro.h
  @author Rahul S. Sampath
  */

#ifndef __DENDRO_H__
#define __DENDRO_H__

#ifdef __USE_64_BIT_INT__
#define DendroIntL long long
#define DendroIntLSpecifier %lld
#define DendroUIntLSpecifier %llu
#else
#define DendroIntL int
#define DendroIntLSpecifier %d
#define DendroUIntLSpecifier %u
#endif

#endif

