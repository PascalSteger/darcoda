module tree_data
 implicit none
 public
 integer,parameter :: bucketmax=1000,maxneigh=64
 real*8,parameter  :: dx0=1.0,dither=-0.0
 character*72::filename
 real*8, parameter::unit_l=  0.261427528685697D+24
 real*8, parameter::mH    =  1.667d-24
 real*8, parameter::msun  =  1.989d33
end module tree_data

module tdef
  implicit none
  type llist
    integer index
    type(llist), pointer::next
  end type llist

  type item
    integer id, level, iparent,ichild
    integer, dimension(33)::neighbor
    real*8 x,y,z,mass,rho,h
  end type item

  type node
    integer index,level
    integer nm0,nm1,nm2,nm3,nm4,nm5,nm6,nm7,nm8
    real*8::down,left,front,xc,yc,zc,dx
    type(node), pointer::c1,c2,c3,c4,c5,c6,c7,c8,parent
    type(llist),pointer::m0,m1,m2,m3,m4,m5,m6,m7,m8
  end type node

  public:: llist,node,item

  contains

recursive subroutine insert_item (index, root)
    implicit none
    integer index
    type(llist), pointer :: root
    if(.not.associated(root)) then
      allocate(root)
      root%index=index
      nullify(root%next)   
    else
      call insert_item(index,root%next)
    endif
end subroutine

subroutine new_node(id,level,parent,down,left,front,dx)
    implicit none
    integer id,level
    type(node), pointer :: parent
    real*8::front,down,left,right,dx
    if(.not.associated(parent))then
      allocate(parent)
      parent%index=id
      parent%level=level
      parent%xc=left+dx*.5
      parent%yc=down+dx*.5
      parent%zc=front+dx*.5
      parent%dx=dx
      parent%down=down
      parent%left=left
      parent%front=front
      parent%nm1=0;parent%nm2=0;parent%nm3=0;parent%nm4=0;parent%nm0=0
      parent%nm5=0;parent%nm6=0;parent%nm7=0;parent%nm8=0
      nullify(parent%parent)
      nullify(parent%c1)
      nullify(parent%c2)
      nullify(parent%c3)
      nullify(parent%c4)
      nullify(parent%c5)
      nullify(parent%c6)
      nullify(parent%c7)
      nullify(parent%c8)
      nullify(parent%m0)
      nullify(parent%m1)
      nullify(parent%m2)
      nullify(parent%m3)
      nullify(parent%m4)
      nullify(parent%m5)
      nullify(parent%m6)
      nullify(parent%m7)
      nullify(parent%m8)
    else
      print *, " Expecting new node, but pointer is already allocated."
      print *, " Stopping."
      stop
    endif
end subroutine

recursive subroutine insert_tree_item(part,parent)
    implicit none
    type(item)::part
    type(node), pointer :: parent
    real*8::up,down,left,front,dx,x,y,z,dxh
    x=part%x
    y=part%y
    z=part%z
    down=parent%down
    left=parent%left
    front=parent%front
    dx=parent%dx
    dxh=dx*.5

    if(.not.associated(parent%m0))then
      allocate(parent%m0)
      parent%nm0=parent%nm0+1
      nullify(parent%m0%next)
    else
      parent%nm0=parent%nm0+1 
      call insert_item(part%id,parent%m0%next)
    endif

    if(x<left+dxh .and. x>left)then
     if (y<down+dxh .and. y > down)then
      if(z<front+dxh .and. z > front) then
       if(.not.associated(parent%m3))then
         allocate(parent%m3)
         parent%m3%index=part%id
         parent%nm3=parent%nm3+1
         nullify(parent%m3%next)
       else
         parent%nm3=parent%nm3+1
         call insert_item(part%id,parent%m3%next)
       endif
      else if(z<front+dx .and. z>front) then
       if(.not.associated(parent%m7))then
         allocate(parent%m7)
         parent%m7%index=part%id
         parent%nm7=parent%nm7+1
         nullify(parent%m7%next)
       else
         parent%nm7=parent%nm7+1
         call insert_item(part%id,parent%m7%next)
       endif
      endif
     else if (y<down+dx .and. y>down) then
      if(z<front+dxh .and. z> front) then
       if(.not.associated(parent%m2))then
         allocate(parent%m2)
         parent%m2%index=part%id
         parent%nm2=parent%nm2+1
         nullify(parent%m2%next)
       else
         parent%nm2=parent%nm2+1
         call insert_item(part%id,parent%m2%next)
       endif
      else if(z<front+dx .and. z>front) then
       if(.not.associated(parent%m6))then
         allocate(parent%m6)
         parent%m6%index=part%id
         parent%nm6=parent%nm6+1
         nullify(parent%m6%next)
       else
         parent%nm6=parent%nm6+1
         call insert_item(part%id,parent%m6%next)
       endif
      endif
     endif
    else if(x<left+dx .and. x>left) then
      if(y<down+dxh .and. y>down)then
       if(z<front+dxh .and. z>front) then
        if(.not.associated(parent%m4))then
         allocate(parent%m4)
         parent%m4%index=part%id
         parent%nm4=parent%nm4+1
         nullify(parent%m4%next)
        else
         parent%nm4=parent%nm4+1
         call insert_item(part%id,parent%m4%next)
        endif
       else if(z<front+dx .and. z>front) then
        if(.not.associated(parent%m8))then
         allocate(parent%m8)
         parent%m8%index=part%id
         parent%nm8=parent%nm8+1
         nullify(parent%m8%next)
        else
         parent%nm8=parent%nm8+1
         call insert_item(part%id,parent%m8%next)
        endif
       endif
      else if (y<down+dx .and. y > down) then
       if(z<front+dxh .and. z>front) then
        if(.not.associated(parent%m1))then
         allocate(parent%m1)
         parent%m1%index=part%id
         parent%nm1=parent%nm1+1
         nullify(parent%m1%next)
        else
         parent%nm1=parent%nm1+1
         call insert_item(part%id,parent%m1%next)
        endif
       else if(z<front+dx .and. z>front) then
        if(.not.associated(parent%m5))then
         allocate(parent%m5)
         parent%m5%index=part%id
         parent%nm5=parent%nm5+1
         nullify(parent%m5%next)
        else
         parent%nm5=parent%nm5+1
         call insert_item(part%id,parent%m5%next)
        endif
       endif
      endif
    endif 
end subroutine

recursive subroutine build_tree(id,level,nmax,npart,part,parent)
   implicit none
   integer ichild,nmax,n,npart,id,level
   real*8 bucket,down,left,front,dx
   type(item),dimension(:)::part
   type(node), pointer::parent,child
   type(llist),pointer::current
 
   dx=parent%dx*.5
   do ichild=1,8
    select case(ichild)
    case(1)
     n=parent%nm1
     down=parent%down+dx
     left=parent%left+dx
     front=parent%front
    case(2)
     n=parent%nm2
     down=parent%down+dx
     left=parent%left
     front=parent%front
    case(3)
     n=parent%nm3
     down=parent%down
     left=parent%left
     front=parent%front
    case(4)
     n=parent%nm4
     down=parent%down
     left=parent%left+dx
     front=parent%front
    case(5)
     n=parent%nm5
     down=parent%down+dx
     left=parent%left+dx
     front=parent%front+dx
    case(6)
     n=parent%nm6
     down=parent%down+dx
     left=parent%left
     front=parent%front+dx
    case(7)
     n=parent%nm7
     down=parent%down
     left=parent%left
     front=parent%front+dx
    case(8)
     n=parent%nm8
     down=parent%down
     left=parent%left+dx
     front=parent%front+dx
    end select
    if (n>nmax) then
      id=id+1
      level=parent%level+1
      select case(ichild)
      case(1)
        nullify(child)
        call new_node(id,level,child,down,left,front,dx)
        child%parent=>parent
        parent%c1=>child
        nullify(current)
        current=>parent%m1
        do while (associated(current))
          call insert_tree_item(part(current%index),child)
          current=>current%next
        enddo
        call build_tree(id,level,nmax,npart,part,child)
      case(2)
        nullify(child)
        call new_node(id,level,child,down,left,front,dx)
        child%parent=>parent
        parent%c2=>child
        nullify(current)
        current=>parent%m2
        do while (associated(current))
          call insert_tree_item(part(current%index),child)
          current=>current%next
        enddo
        call build_tree(id,level,nmax,npart,part,child)
      case(3)
        nullify(child)
        call new_node(id,level,child,down,left,front,dx)
        child%parent=>parent
        parent%c3=>child
        nullify(current)
        current=>parent%m3
        do while (associated(current))
          call insert_tree_item(part(current%index),child)
          current=>current%next
        enddo
        call build_tree(id,level,nmax,npart,part,child)
      case(4)
        nullify(child)
        call new_node(id,level,child,down,left,front,dx)
        child%parent=>parent
        parent%c4=>child
        nullify(current)
        current=>parent%m4
        do while (associated(current))
          call insert_tree_item(part(current%index),child)
          current=>current%next
        enddo
        call build_tree(id,level,nmax,npart,part,child)
      case(5)
        nullify(child)
        call new_node(id,level,child,down,left,front,dx)
        child%parent=>parent
        parent%c5=>child
        nullify(current)
        current=>parent%m5
        do while (associated(current))
          call insert_tree_item(part(current%index),child)
          current=>current%next
        enddo
        call build_tree(id,level,nmax,npart,part,child)
      case(6)
        nullify(child)
        call new_node(id,level,child,down,left,front,dx)
        child%parent=>parent
        parent%c6=>child
        nullify(current)
        current=>parent%m6
        do while (associated(current))
          call insert_tree_item(part(current%index),child)
          current=>current%next
        enddo
        call build_tree(id,level,nmax,npart,part,child)
      case(7)
        nullify(child)
        call new_node(id,level,child,down,left,front,dx)
        child%parent=>parent
        parent%c7=>child
        nullify(current)
        current=>parent%m7
        do while (associated(current))
          call insert_tree_item(part(current%index),child)
          current=>current%next
        enddo
        call build_tree(id,level,nmax,npart,part,child)
      case(8)
        nullify(child)
        call new_node(id,level,child,down,left,front,dx)
        child%parent=>parent
        parent%c8=>child
        nullify(current)
        current=>parent%m8
        do while (associated(current))
          call insert_tree_item(part(current%index),child)
          current=>current%next
        enddo
        call build_tree(id,level,nmax,npart,part,child)
      end select
    endif
  enddo

end subroutine

recursive subroutine id_tree(nid,maxlevel,npart,part,parent)
   implicit none
   integer ichild,n,npart,nid,maxlevel
   type(item),dimension(:):: part
   type(node), pointer::parent,child
   type(llist), pointer::current,member

   nullify(child)
   nullify(member)
   
   do ichild=1,8
    select case(ichild)
    case(1)
      child=>parent%c1
      member=>parent%m1
      n=parent%nm1
    case(2)
      child=>parent%c2
      member=>parent%m2
      n=parent%nm2
    case(3)
      child=>parent%c3
      member=>parent%m3
      n=parent%nm3
    case(4)
      child=>parent%c4
      member=>parent%m4
      n=parent%nm4
    case(5)
      child=>parent%c5
      member=>parent%m5
      n=parent%nm5
    case(6)
      child=>parent%c6
      member=>parent%m6
      n=parent%nm6
    case(7)
      child=>parent%c7
      member=>parent%m7
      n=parent%nm7
    case(8)
      child=>parent%c8
      member=>parent%m8
      n=parent%nm8
    end select
    if(associated(child))then
     call id_tree(nid,maxlevel,npart,part,child)
    else
     if(n>0)then
       nullify(current)
       current=>member
       do while(associated(current))
         maxlevel=max(maxlevel,parent%level)
         part(current%index)%level=parent%level
         part(current%index)%iparent=parent%index
         part(current%index)%ichild=ichild
         current=>current%next
       enddo
     endif
    endif
   enddo     
end subroutine

recursive subroutine destroy_tree(parent)
   implicit none
   integer ichild,matchId
   type(node), pointer::parent,child
   logical found_match

   nullify(child)
   
   do ichild=1,8
    select case(ichild)
    case(1)
      child=>parent%c1
    case(2)
      child=>parent%c2
    case(3)
      child=>parent%c3
    case(4)
      child=>parent%c4
    case(5)
      child=>parent%c5
    case(6)
      child=>parent%c6
    case(7)
      child=>parent%c7
    case(8)
      child=>parent%c8
    end select
    if(associated(child))then
       call destroy_tree(child)
    else
       deallocate(parent)
       exit
    endif
   enddo     
end subroutine

recursive subroutine print_tree(npart,part,parent)
   implicit none
   integer ichild,n,npart
   type(item),dimension(:):: part
   type(node), pointer::parent,child
   type(llist), pointer::current,member

   nullify(child)
   nullify(member)
   
   do ichild=1,8
    select case(ichild)
    case(1)
      child=>parent%c1
      member=>parent%m1
      n=parent%nm1
    case(2)
      child=>parent%c2
      member=>parent%m2
      n=parent%nm2
    case(3)
      child=>parent%c3
      member=>parent%m3
      n=parent%nm3
    case(4)
      child=>parent%c4
      member=>parent%m4
      n=parent%nm4
    case(5)
      child=>parent%c5
      member=>parent%m5
      n=parent%nm5
    case(6)
      child=>parent%c6
      member=>parent%m6
      n=parent%nm6
    case(7)
      child=>parent%c7
      member=>parent%m7
      n=parent%nm7
    case(8)
      child=>parent%c8
      member=>parent%m8
      n=parent%nm8
    end select
    if(associated(child))then
!     print *, child%xc,child%yc,child%zc, &
!              ichild,child%index,child%level,n
     call print_tree(npart,part,child)
    else
     if(n>0)then
       nullify(current)
       current=>member
       do while(associated(current))
         print *, part(current%index)%x,part(current%index)%y,&
                  part(current%index)%z,   &
                  part(current%index)%mass,&
                  part(current%index)%rho, &
                  part(current%index)%h,part(current%index)%level,&
                  part(current%index)%iparent,part(current%index)%ichild
         current=>current%next
       enddo
       print *, "#Done with current leaf"
     endif
    endif
   enddo     
end subroutine

recursive subroutine get_neighbors(maxlevel,maxneigh,npart,part,parent)
   implicit none
   integer npart,i,maxneigh,maxlevel,ichild,thischild,nlocal
   type(item),dimension(:):: part
   type(node), pointer::parent,child
   type(llist), pointer::current,member

   print *,'#in get_neighbors'
   nullify(child)
   nullify(member)

   if(parent%level==maxlevel)then
      print *,'#parentlevel=maxlevel'
      thischild=0
      member=>parent%m0
      print *,'#member link set'
      nlocal=parent%nm0
      print *,'#nlocal set: ',nlocal
      !print *,"#", ichild,child%level,maxlevel,parent%level,parent%index
      !print *,"#", parent%level,maxlevel,parent%index
      !print *,"# located ",parent%xc,parent%yc,parent%zc,parent%dx
      print *,'#calling gather'
      call gather(thischild,maxneigh,nlocal,npart,part,member,parent)
   elseif(parent%level<maxlevel)then
      print *,'# step through 8 childs'
      do ichild=1,8
         select case(ichild)
         case(1)
            child=>parent%c1
            member=>parent%m1
            nlocal=parent%nm1
         case(2)
            child=>parent%c2
            member=>parent%m2
            nlocal=parent%nm2
         case(3)
            child=>parent%c3
            member=>parent%m3
            nlocal=parent%nm3
         case(4)
            child=>parent%c4
            member=>parent%m4
            nlocal=parent%nm4
         case(5)
            child=>parent%c5
            member=>parent%m5
            nlocal=parent%nm5
         case(6)
            child=>parent%c6
            member=>parent%m6
            nlocal=parent%nm6
         case(7)
            child=>parent%c7
            member=>parent%m7
            nlocal=parent%nm7
         case(8)
            child=>parent%c8
            member=>parent%m8
            nlocal=parent%nm8
         end select
         print *,'#  end select'
         if(associated(child))then
            ! guess: never called, from elseif in upper level
            !if(parent%level>maxlevel)then
            !   print *,'#error! should never get called'
            !   thischild=0
            !   member=>parent%m0
            !   nlocal=parent%nm0
            !   print *, ichild,child%level,maxlevel,parent%level,parent%index
            !   print *," located ",parent%xc,parent%yc,parent%zc,parent%dx
            !   call gather(thischild,maxneigh,npart,part,member,parent)
            !else
               call get_neighbors(maxlevel,maxneigh,npart,part,child)
            !endif
         else
            print *,'# found child without children'
            ! found child without children
            ! try to calculate neighbors for all particles 
            ! associated with this child
            ! only descend to maxlevel
            call gather(ichild,maxneigh,nlocal,npart,part,member,parent)
         endif
      enddo
   endif
 end subroutine get_neighbors

recursive subroutine gather(ichild,maxneigh,nlocal,npart,part,member,parent)
   implicit none
   integer npart,maxneigh,ineigh,i_out,i_in,ichild,dist2maxId,i,nlocal
   real*8 dist2max,xo,yo,zo,xi,yi,zi,h,dist2,xc,yc,zc,dx,spacing,kernel
   real*8,dimension(maxneigh)::dist2arr,massarr
   integer, dimension(maxneigh)::dist2int
   integer, dimension(:), allocatable:: indexarr
   type(item),dimension(:)::part
   type(item),dimension(:),allocatable::m0part
   type(llist), pointer::inner,outer,member
   type(node), pointer::parent
   logical findEdge

   print *,'# in gather'
   xc=parent%xc;yc=parent%yc;zc=parent%zc;dx=parent%dx

   allocate(m0part(1:nlocal))
   allocate(indexarr(1:nlocal))
   nullify(outer)
   outer=>member
   i=0
   do while(associated(outer))
    i=i+1
    i_out=outer%index
    print *, i,nlocal,i_out
    m0part(i)=part(i_out)
    indexarr(i)=i_out
    outer=>outer%next
   enddo
   do i_out=1,nlocal
     h=part(i_out)%h
     if(h<0.)then
       xo=m0part(i_out)%x;yo=m0part(i_out)%y;zo=m0part(i_out)%z
       dist2int=0
       dist2arr=0.

       ineigh=1
       do i_in=1,nlocal
         if(i_out/=i_in)then
            xi=m0part(i_in)%x
            yi=m0part(i_in)%y
            zi=m0part(i_in)%z
            dist2=(xo-xi)**2+(yo-yi)**2+(zo-zi)**2
            if(ineigh<maxneigh)then
              massarr(ineigh)=part(i_in)%mass
              dist2arr(ineigh)=dist2
              dist2int(ineigh)=indexarr(i_in)
              ineigh=ineigh+1
            elseif(ineigh==maxneigh)then
              dist2arr(ineigh)=dist2
              dist2int=i_in
              massarr(ineigh)=part(i_in)%mass
              dist2max=0.
              dist2maxId=1
              do i=1,ineigh
                if(dist2arr(i)>dist2max)then
                  dist2max=dist2arr(i)
                  dist2maxId=i
                endif
              enddo
              ineigh=ineigh+1
            else
              if (dist2<dist2max)then
                dist2arr(dist2maxId)=dist2
                dist2int(dist2maxId)=i_in
                dist2max=0.
                dist2maxId=1
                i=1
                do while (i<ineigh)
                  if(dist2arr(i)>dist2max)then
                    dist2max=dist2arr(i)
                    dist2maxId=i
                  endif
                  i=i+1
                enddo
              endif
            endif
          endif
        enddo
	h=sqrt(dist2max)
        if(.not.findEdge(ichild,xo,yo,zo,h,xc,yc,zc,dx,parent%level).and. &
        ineigh>maxneigh)then
	  part(indexarr(i_out))%h=h
          part(indexarr(i_out))%rho=kernel(maxneigh,h,massarr,dist2arr)
          part(indexarr(i_out))%neighbor=dist2int
        endif
     endif ! smoothing length if statement
   enddo
   deallocate(indexarr,m0part)
end subroutine

subroutine top_gather(maxneigh,ipart,npart,part)
    implicit none 
    integer ipart,i,innerpart
    integer nneigh,maxneigh,dist2maxId,ichild,npart
    real*8 dist2max,xo,yo,zo,xi,yi,zi,h,dist2,xc,yc,zc,dx,kernel
    real*8,dimension(maxneigh)::dist2arr,dist2int,massarr
    type(item),dimension(:)::part

    h=part(ipart)%h
      xo=part(ipart)%x
      yo=part(ipart)%y
      zo=part(ipart)%z
        nneigh=1.
        dist2int=0
        dist2arr=0.

        do innerpart=1,npart
          if(innerpart/=ipart)then
            xi=part(innerpart)%x
            yi=part(innerpart)%y
            zi=part(innerpart)%z
            dist2=(xo-xi)**2+(yo-yi)**2+(zo-zi)**2
            if(nneigh<maxneigh)then
              massarr(nneigh)=part(innerpart)%mass
              dist2arr(nneigh)=dist2
              dist2int(nneigh)=innerpart
              nneigh=nneigh+1
            elseif(nneigh==maxneigh)then
              dist2arr(nneigh)=dist2
              dist2int=innerpart
              massarr(nneigh)=part(innerpart)%mass
              dist2max=0.
              dist2maxId=1
              do i=1,nneigh
                if(dist2arr(i)>dist2max)then
                  dist2max=dist2arr(i)
                  dist2maxId=i
                endif
              enddo
              nneigh=nneigh+1
            else
              if (dist2<dist2max)then
                dist2arr(dist2maxId)=dist2
                dist2int(dist2maxId)=innerpart
                dist2max=0.
                dist2maxId=1
                i=1
                do while (i<nneigh)
                  if(dist2arr(i)>dist2max)then
                    dist2max=dist2arr(i)
                    dist2maxId=i
                  endif
                  i=i+1
                enddo
              endif
            endif
          endif
        enddo
	h=sqrt(dist2max)
        part(ipart)%h=h
        part(ipart)%rho=kernel(maxneigh,h,massarr,dist2arr)
        part(ipart)%neighbor=dist2int
        print *,'#end of gather'
end subroutine
end module tdef


program main
 use tdef
 use tree_data
 implicit none
 integer :: i,id, level,thisneigh,neigh,id1,level1,maxlevel,ilevel,&
            ichild,idither,ierror,iheader,nentry
 real*8 :: val
 real*8 :: bucket,newtime,oldtime,left,down,front
 type(item), dimension(:), allocatable::point
 type(node),  pointer::tree
 character*2 ch

 nullify(tree)

 iheader=-2 ! last two entries should not be counted
 nentry=0
 call getarg(1,filename)
 print *, '#',filename
 open(unit=11,file=filename,form="formatted")
 do 
  read(11,"(2A)",iostat=ierror)ch
  if(ierror<0)exit
  if(ch(1:1)=='#'.or.ch(2:2)=='#')then
    iheader=iheader+1
  else
    nentry=nentry+1
  endif
 enddo
 if (ierror>0)then
   print *, "Error occurred during read"
   stop
 endif
 if (nentry==0)then
   print *, "No valid record line"
   stop
 endif 

 print *, "# allocating memory for N = ",nentry," particles."
 allocate(point(1:nentry)) 

 rewind(11)
 do i=1,iheader
  read(11,*)
 enddo
 do i=1,nentry
  read(11,*) point(i)%mass, point(i)%x,point(i)%y,point(i)%z
!  print *, point(i)%mass, point(i)%x,point(i)%y, point(i)%z
  point(i)%id=i
 enddo
 close(11)
 
 do i=1,nentry
   point(i)%h=-1.
 enddo
 print *, "#Set smoothing length"

 left=-dither;down=-dither;front=-dither

 call cpu_time(oldtime)
 id=1; level=1
 call new_node(id,level,tree,down,left,front,dx0) ! minus gives root node id=0
 do i=1,nentry
   call insert_tree_item(point(i),tree)
 enddo  
 if(associated(tree%c1))tree%c1%parent=>tree
 if(associated(tree%c2))tree%c2%parent=>tree
 if(associated(tree%c3))tree%c3%parent=>tree
 if(associated(tree%c4))tree%c4%parent=>tree
 if(associated(tree%c5))tree%c5%parent=>tree
 if(associated(tree%c6))tree%c6%parent=>tree
 if(associated(tree%c7))tree%c7%parent=>tree
 if(associated(tree%c8))tree%c8%parent=>tree
 
 print *, "# Created root node"

 ! we now have a top-level node with four children.  
 ! We now need to test whether
 ! each child needs to become a parent node, too.
 
 call build_tree(id,level,bucketmax,nentry,point,tree)
 call cpu_time(newtime)
 print *,"#Built main tree",newtime-oldtime
 oldtime=newtime

 maxlevel=1
 call id_tree(id,maxlevel,nentry,point,tree)
 call cpu_time(newtime)
 print *, "#id tree found maxlevel ",maxlevel,newtime-oldtime
 oldtime=newtime
 !maxlevel=maxlevel-1
 do ilevel=maxlevel,2,-1
   call get_neighbors(ilevel,maxneigh,nentry,point,tree)
 enddo
 call cpu_time(newtime)
 print *, "# Through tree except trunk ",newtime-oldtime
 oldtime=newtime
 do i=1,nentry
  if(point(i)%h<0.) call top_gather(maxneigh,i,nentry,point)
 enddo
 call cpu_time(newtime)
 print *, "# Done with trunk ",newtime-oldtime

 do i=1,nentry ! convert to #/cc
   point(i)%rho=point(i)%rho*msun/(mH*unit_l**3)
 enddo

 call print_tree(nentry,point,tree)
 print *,'#tree printed'
 call destroy_tree(tree)
 print *,'#tree destroyed'
 deallocate(point)
 print *,'#point deallocated'

end program
 
logical function findEdge(nchild,x,y,z,h,xc,yc,zc,dx,level)
    implicit none
    integer nchild,level
    real*8 x,y,z,h,l,f,d,dxh
    real*8, intent(in)::xc,yc,zc,dx

     findEdge=.false.
     if(level==1)return
     dxh=dx*0.5
     select case(nchild)
     case(0)
       l=xc-dxh
       d=yc-dxh
       f=zc-dxh
       dxh=dx
     case(1)
       l=xc
       d=yc
       f=zc-dxh
     case(2)
       l=xc-dxh
       d=yc
       f=zc-dxh
     case(3)
       l=xc-dxh
       d=yc-dxh
       f=zc-dxh
     case(4)
       l=xc
       d=yc-dxh
       f=zc-dxh
     case(5)
       l=xc
       d=yc
       f=zc
     case(6)
       l=xc-dxh
       d=yc
       f=zc
     case(7)
       l=xc-dxh
       d=yc-dxh
       f=zc
     case(8)
       l=xc
       d=yc-dxh
       f=zc
     end select

     if(2.*h>abs(x-l      ))findEdge=.true.
     if(2.*h>abs(x-(l+dxh)))findEdge=.true.
     if(2.*h>abs(y-d      ))findEdge=.true.
     if(2.*h>abs(y-(d+dxh)))findEdge=.true.
     if(2.*h>abs(z-f      ))findEdge=.true.
     if(2.*h>abs(z-(f+dxh)))findEdge=.true.

!     print *, level,nchild,down,left,d,l
    return
end function

real*8 function kernel(n,h,mass,dist)
    implicit none
    integer i,n
    real*8 q,h
    real*8,dimension(n)::mass,dist

    dist=sqrt(dist)
    kernel=0.

    do i=1,n
      q=dist(i)/h
      if(q<1.)then
        kernel=kernel+(1.+q*(q*(-1.5+.75*q)))*mass(i)
      else if(q<2.)then
        kernel=kernel+.25*(2.-q)**3*mass(i)
      endif
    enddo
    kernel=kernel/(h**3*3.1415926535897931d0)
    return
end function

