#!/usr/bin/env python3

##
# @file
# test exception handling

def fa():
    raise Exception('negative','found')
    print('this should not be printed!')
    return 3

def fb():
    return fa()

def fc():
    return fb()

print('if no line is printed after this, the program has exited too fast..')
try:
    fc()
except Exception as detail:
    print('Handling run-time error:', detail)

print('no! exception handled,  execution has beed resumed')
