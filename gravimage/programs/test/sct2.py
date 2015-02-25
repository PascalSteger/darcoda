from mpi4py import MPI

comm = MPI.COMM_WORLD

rank = comm.rank
size = comm.size

print('Hi, my rank is:', rank)

if comm.rank == 0:
    print("doing the task of rank 0")


if comm.rank == 1:
    print("doing the task of rank 1")
    print("rank:"+str(rank))
    print("size:"+str(size))
    print(9**(rank+3))

if comm.rank == 2:
    print("doing the task of rank 2")
    print("rank:"+str(rank))
    print("size:"+str(size))
    print(9**(rank+3))
