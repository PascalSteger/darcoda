#!/usr/bin/python

import numpy
import array
import struct
import sys

VERBOSE=1
URM=4

def getFortranFile(file,s):
	a=numpy.fromfile(file,numpy.dtype(s))
        return a

def getFortranData(N,RM,file,doArray,type):
        if RM == 4:
                SRM = "i"
        elif RM == 8:
                SRM = "l"
        else:
                print "Your choice for RM is unreasonable."
                print "Record marker should be 4 or 8."
                sys.exit(1)
        if type == "int":
                t='i'
		RT=4
        elif type == "float":
                t='f'
		RT=4
	elif type == "double":
		t='d'
		RT=8
        else:
                print "Your choice for type is unreasonable."
                print "Type is limited to float or int."
                print "Float works for doubles, too."
                sys.exit(1)
        if doArray==0:
                s=file.read(RM)
		s=file.read(RT)
		val, = struct.unpack(t,s)
		s=file.read(RM)
        elif doArray==1:
                s = file.read(RM)
		RL,=struct.unpack(SRM,s)
		if RL/RT != N:
			print "Record length and array size are not congruent.", RL/RT,N
			print "Maybe dealing with really large arrays?"
                binvalues=array.array(t)
                binvalues.read(file,N)
                val=numpy.array(binvalues,type)
		s = file.read(RM)
        else:
                print "Your choice for doArray is unreasonable."
                print "Print use 0 for single value and 1 for array."
                sys.exit(1)
                
        return val

def getJunkRead(RM,file):
        if RM == 4:
                SRM = "i"
        elif RM == 8:
                SRM = "l"
        else:
                print "Your choice for RM is unreasonable."
                print "Record marker should be 4 or 8."
                sys.exit(1)
        s=file.read(RM)
	var,=struct.unpack(SRM,s)
        s=file.read(var)
        s=file.read(RM)
	return var

def getRamsesPartHead(noutput,dir):
	"""Return ncpu, ndim, and npart from fortran RAMSES data"""
	soutput=repr(noutput).zfill(5)
	dir+=soutput
	filename=dir+"/part/part_"+soutput+".out00001"
	if VERBOSE==1:print "#",filename
	try:
		f=open(filename,"rb")
	except IOError:
		print "Problem with file"+filename+". Shutting down."
		sys.exit(1)
	ncpu=getFortranData(1,URM,f,0,"int")
	ndim=getFortranData(1,URM,f,0,"int")
	f.close()
		
	npart=0
 	for i in xrange(ncpu):
		scpu=repr(i+1).zfill(5)
		filename=dir+"/part/part_"+soutput+".out"+scpu
		if VERBOSE==1: print "#",filename
		f=open(filename,"rb")
		idum=getFortranData(1,URM,f,0,"int")
		idum=getFortranData(1,URM,f,0,"int")
		idum=getFortranData(1,URM,f,0,"int")
		npart+=idum
		f.close()
	return ncpu,ndim,npart

def getRamsesInfo(noutput,dir):
        """Return boxlen,time,h0,aexp,omega_m,omega_l,omega_k,omega_b,unit_l,unit_d,unit_t from RAMSES data"""
        soutput=repr(noutput).zfill(5)
        dir+=soutput
	filename=dir+"/info_"+soutput+".txt"
        if VERBOSE==1:print "#",filename
        try:
                f=open(filename,"rb")
        except: # [PS] except IOError:
                print "Problem with file"+f+". Shutting down."
                sys.exit(1)

	for line in f:
		try:
			key,val=line.split("=")
		except:
			key=" ";val=" "
		key=key.strip()
		val=val.strip()
		if key=="ncpu": ncpu=int(val)
		elif key=="ndim": ndim=int(val)
		elif key=="levelmin": levelmin=int(val)
		elif key=="levelmax": levelmax=int(val)
		elif key=="ngridmax": ngridmax=int(val)
		elif key=="nstep_coarse": nstep_coarse=int(val)
		elif key=="boxlen": boxlen=float(val)
		elif key=="time": time=float(val)
		elif key=="aexp": aexp=float(val)
		elif key=="H0": H0=float(val)
		elif key=="omega_m": omega_m=float(val)
		elif key=="omega_l": omega_l=float(val)
		elif key=="omega_k": omega_k=float(val)
		elif key=="omega_b": omega_b=float(val)
		elif key=="unit_l": unit_l=float(val)
		elif key=="unit_d": unit_d=float(val)
		elif key=="unit_t": unit_t=float(val)
	f.close()
	return ncpu,ndim,levelmin,levelmax,ngridmax,nstep_coarse,boxlen,time,aexp,H0,omega_m,omega_l,omega_k,omega_b,unit_l,unit_d,unit_t

def getRamsesPart(npart,ndim,ncpu,noutput,dir):
	"""Read RAMSES particle data. Returns X, V, Mass, and ID"""
        soutput=repr(noutput).zfill(5)
        dir+=soutput
	x=numpy.empty((ndim,npart),'float')
	v=numpy.empty((ndim,npart),'float')
	mass=numpy.empty((npart),'float')
	id=numpy.empty((npart),'int')
	
	ii=0
	tpart=0
        for i in xrange(ncpu):
                scpu=repr(i+1).zfill(5)
                filename=dir+"/part/part_"+soutput+".out"+scpu
                if VERBOSE==1: print "#",filename
                f=open(filename,"rb")
		ncpu=getFortranData(1,URM,f,0,"int")
		ndim=getFortranData(1,URM,f,0,"int")
		ipart=getFortranData(1,URM,f,0,"int")
		idum=getJunkRead(4,f)
                idum=getJunkRead(4,f)
                idum=getJunkRead(4,f)
                idum=getJunkRead(4,f)
                idum=getJunkRead(4,f)
		tmpX=numpy.empty((ndim,ipart),"float")
		tmpV=numpy.empty((ndim,ipart),"float")
		for idum in xrange(ndim):
			tmpX[idum][:]=getFortranData(ipart,URM,f,1,"double")
		for idum in xrange(ndim):
			tmpV[idum][:]  =getFortranData(ipart,URM,f,1,"double")
		tmpMass=getFortranData(ipart,URM,f,1,"double")
		tmpId  =getFortranData(ipart,URM,f,1,"int")
		for j in xrange(ipart):
			for idum in xrange(ndim):
				x[idum][ii]=tmpX[idum][j]	
				v[idum][ii]=tmpV[idum][j]	
			mass[ii]=tmpMass[j]
			id[ii]=tmpId[j]
			ii+=1
		f.close()
			
	return x,v,mass,id

def getRamsesKey(noutput,dir):
        """Return boxlen,time,h0,aexp,omega_m,omega_l,omega_k,omega_b,unit_l,unit_d,unit_t from RAMSES data"""
        soutput=repr(noutput).zfill(5)
        dir+=soutput
	filename=dir+"/info_"+soutput+".txt"
        if VERBOSE==1:print "#",filename
        try:
                f=open(filename,"rb")
        except IOerror:
                print "Problem with file"+f+". Shutting down."
                sys.exit(1)
	lines=f.readlines()
	nLines=len(lines)
	for i in xrange(nLines):
		try:
			key,val=line[i].split("=")
		except:
			key=" ";val=" "
		key=key.strip()
		val=val.strip()
		if key=="ncpu": ncpu=int(val)
		if lines.find("hilbert")>-1:
			i+=2
			for icpu in xrange(ncpu):
				icheck,bound_key[icpu],bound_key[icpu+1]=lines[i].split()
				bound_key=numpy.empty(ncpu+1,'int')
				cpu_read=numpy.empty(ncpu,'int')
				if icheck!=icpu:
					print "Hilbert key and cpus do not match."
					sys.exit(1)
		

def swap(a,b):
      
      swapit = a
      a = b
      b = swapit
      
      return a,b


def shiftdown(a,ia,start,end,n):

      root = start

      while (root*2+1 <= end):
            child = root*2+1
            if (child < end and a[child] < a[child+1]): child=child+1
            if (a[root] < a[child]):
                  a[root],a[child]=swap(a[root],a[child])
                  ia[root],ia[child]=swap(ia[root],ia[child])
            root = child

      return a,ia

def heapit(a,ia,count,n):

      start = count/2-1

      while start >= 0:
            a,ia=shiftdown(a,ia,start,count-1,n)
            start = start-1
      return a,ia

def heapsort(b,ib,count,n):
      a=b
      ia=ib

      a,ia=heapit(a,ia,count,n)

      end = count-1
      while end > 0:
            a[end],a[0]=swap(a[end],a[0])
            ia[end],ia[0]=swap(ia[end],ia[0])
            end = end -1
            a,ia=shiftdown(a,ia,0,end,n)
      return a,ia
