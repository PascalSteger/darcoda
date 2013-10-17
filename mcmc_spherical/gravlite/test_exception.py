#!/usr/bin/python3.3

def fa():
    raise Exception('negative','found')
    print('this should not be printed!')
    return 3

def fb():
    return fa()

def fc():
    return fb()

try:
    fc()
except Exception as detail:
    print('Handling run-time error:', detail)

print('if this line is printed, execution was resumed after handling')
