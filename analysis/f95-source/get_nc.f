       subroutine get_nc(nres,imap,kcont,klist,ch,inn)
       implicit none
       integer nres,kcont,inn
       integer imap(nres,nres)
       character*1 ch(nres)
       integer klist(kcont,2)
       integer i,j,k
       
       inn=0
       do k=1,kcont
        i=klist(k,1)
        j=klist(k,2)
        if (ch(i).ne.ch(j)) then
         if (imap(i,j).gt.0) then
          inn=inn+1
         endif
        endif
       enddo
       return
       end
       
c
c      subroutine get_nc(nres,imap,nch,ichain,inn)
c      implicit none
c      integer nres
c      integer imap(nres,nres)
c      integer ichain(nch,2)
c      integer nch,inn,inc
c      integer k,i1,i2,j1,j2,j
c      
c      inc=0
c      do i1=1,nres
c       do i2=1,nres
c        if (imap(i1,i2).gt.0) then
c         inc=inc+1
c        endif
c       enddo
c      enddo
c      
c      j=0
c      do k=1,nch
c       i1=ichain(k,1)
c       i2=ichain(k,2)
c       do j1=i1,i2
c        do j2=i1,i2
c         if (imap(j1,j2).gt.0) then
c          j=j+1
c         endif
c        enddo
c       enddo
c      enddo
c      inn=inc-j
c      return
c      end
c