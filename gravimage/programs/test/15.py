#!/usr/bin/env ipython2

def mpi_run():
    if MPI.COMM_WORLD.Get_rank() == 0:
        mpi_controller()
    else:
        mpi_worker()

def mpi_controller():
    ''' 
    Controls the distribution of data-sets to the nodes
    '''
    
    iterations = 10000000
    burnin = 4000000

    orders = range(3, 6)

    # Stores the original task list
    task_list = []
    
    # Stores a list of stats
    stats_list = []
    
    # Stores a list of processes for exit checking
    process_list = range(1, MPI.COMM_WORLD.Get_size())
    
    # Generate the task_list    
    for order in orders:
        for plambda in plambdas:
            task_list.append({'order':order, 'plambda':plambda,
                              'iterations':iterations, 'burnin':burnin})

    total_time_start = time.time()
    print 'Running sampler with: %i tasks on %i processes ' % ( len(task_list), MPI.COMM_WORLD.Get_size())
    
    while len(process_list) > 0:
            
        status = MPI.Status()

        data = MPI.COMM_WORLD.Recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)

        if status.tag > 9:
            if status.tag == 15:
                # Record the data
                stats_list.append(data[0])
        
            task = []
            if len(task_list) > 0:
                task = task_list.pop()
            
            MPI.COMM_WORLD.Send(task, dest=status.source)        

        elif status.tag == 5: # Status
            update_gui()
            pass
        elif status.tag == 2: # Exit
            process_list.remove(status.source)
            print 'Process %i exited' % status.source
        else:
            print 'Unkown tag %i with msg %s' % (status.tag, str(data))

    print "Data Run Finished"
    print "Total Elapsed time: " + str(datetime.timedelta(seconds = (time.time() - total_time_start)))

def mpi_worker():
    '''
    Worker process
    '''
    
    rank = MPI.COMM_WORLD.Get_rank()
    proc_name = MPI.Get_processor_name()

    # Send ready
    MPI.COMM_WORLD.Send([{'rank':rank, 'name':proc_name}], dest=0, tag=10)

    # Start main data loop
    while True:
        # Get some data
        data = MPI.COMM_WORLD.Recv(source=0)
    
        if len(data) == 0: break;

        # Crunch
        time_start = time.time()

        crunch()
        
        time_end = time.time()
        
        # Save results
        results = {'rank':rank, 'name':proc_name, 'start_time': time_start, 'end_time': time_end}
        # Return results
        MPI.COMM_WORLD.Send([results], dest=0, tag=15)
    
    MPI.COMM_WORLD.Send([], dest=0, tag=2)
