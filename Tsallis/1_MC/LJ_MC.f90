module param
implicit none
integer,parameter :: seed=700001 !seed for random numbers

integer,parameter :: ID=1

integer,parameter :: dime=3 !dimension
integer,parameter :: NN=108 !the number of particles
integer,parameter :: MC_max=NN*100000 !upper MC step
integer,parameter :: modRDF=10000000
integer,parameter :: N_mod=NN*10 !output interval
integer,parameter :: outputmod=NN*10000
double precision,parameter :: eta=0.09d0 !displacement coefficient 
double precision,parameter :: length=4.762d0!3.1748d0 !size of box
double precision,parameter :: T=0.3d0 !temperature ,dimensionless


double precision,parameter :: E_hist_min=-900.0d0
double precision,parameter :: E_hist_max=100.0d0
integer,parameter :: Nhist=1000




double precision,parameter :: pi=3.141592653589793d0
integer,parameter :: N_hist=1800
double precision,parameter :: dr=2.0d0*sqrt(3.0d0)*length/dble(N_hist)

end module param
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program LJ_fluid
use param
implicit none
double precision,dimension(dime,NN) :: x !position
double precision,dimension(dime) :: xm !memory x
double precision :: r !distance between particles
double precision :: E_pair !LJ_energy of one pair
double precision :: E_pre !total energy before
double precision :: E_new !total energy after
double precision :: dE !delta total energy
double precision :: dE_hist
integer,dimension(Nhist) :: histogram
integer :: IX !for random
double precision :: RND !random number
integer :: MC_step
integer :: N_select !selected particle number
integer :: i,j,cnt,freeint
integer,dimension(dime) :: l_bc !for boundary condition
integer :: N_A !the number of accepting
integer,dimension(N_hist) :: CNT_RDF
double precision :: coefficient,RRR
character(len=32) :: file1,file2

!open file
write(file1,'("./data",i0,"/Energy.dat")')ID
open(38,file=file1,status='replace')

!insert seed
IX=seed

call outputDETAIL

!initial setting
xm=0.0d0
cnt=0
CNT_RDF=0
N_A=MC_max
histogram=0
dE_hist=(E_hist_max-E_hist_min)/dble(Nhist)


!initial condition
call initial(x)
open(39,file="./initialfcc.dat",status='old')
do i=1,NN
   read(39,*) x(1,i),x(2,i),x(3,i)
end do
close(39)



!energy calc
E_pre=0.0d0
E_new=0.0d0
do i=1,NN
   do j=i+1,NN
      call distance(x,i,j,r)
      call energy(r,E_pair)
      E_new=E_new+E_pair
   end do      
end do

call output(0,x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!        main loop            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do MC_step=1,MC_max
write(*,*) "step=", MC_step

!energy memory
   E_pre=E_new
   E_new=0.0d0

!select particle
   N_select=0
   call random(IX,RND)
   N_select=int(NN*RND)+1

!small displacement
   do i=1,dime
      call random(IX,RND)
      xm(i)=x(i,N_select)
      x(i,N_select)=x(i,N_select)+eta*(2.0d0*RND-1.0d0)
   end do

!boundary condition
   do i=1,dime
      l_bc(i)=int((2.0d0*x(i,N_select)-length)/length)
      x(i,N_select)=x(i,N_select)-length*l_bc(i)
   end do


!energy calc
   E_new=0.0d0
   do i=1,NN
      do j=i+1,NN
         call distance(x,i,j,r)
         call energy(r,E_pair)
         E_new=E_new+E_pair
      end do      
   end do

!delta energy
   dE=0.0d0
   dE=E_new-E_pre

!   write(*,*) "dE,exp=",dE,exp(-dE/T)
write(*,*) "E_new",E_new

!metropolis method
   if(dE>0.0d0)then
      call random(IX,RND)
      if(RND>exp(-dE/T))then
         N_A=N_A-1
         do i=1,dime
            x(i,N_select)=xm(i)
         end do
         E_new=E_pre
         cnt=cnt+1
      end if
   end if          

!output file
   if(mod(MC_step,outputmod)==0)then
      call output(MC_step,x)
   end if
!      
!   if(mod(MC_step,modRDF)==0)then
!      call RDF_calc(x,MC_step,CNT_RDF)
!   end if
!
!histogram
   if(mod(MC_step,1)==0)then
      freeint=int((E_new-E_hist_min)/dE_hist)+1
      if(freeint<=Nhist)then
      histogram(freeint)=histogram(freeint)+1
      end if
   end if



   if(mod(MC_step,N_mod)==0 .and. cnt==0)then
   write(38,*) MC_step/N_mod,E_new
   end if

   cnt=0
end do

!histogram output
write(file2,'("./data",i0,"/hist.dat")')ID
open(41,file=file2,status="replace")
do i=1,Nhist
   write(41,*) dE_hist*(i-0.5d0)+E_hist_min,histogram(i)
end do
close(41)




!total RDF calculation
!open(26,file="./RDF/RDFliquid.dat",status='replace')
!coefficient=4.0*pi*dble(NN)*dr/(length**3.0)*int(MC_max/modRDF)
!do i=1,N_hist
!   RRR=dr*(dble(i)-1.0d0)
!   write(26,*) RRR,CNT_RDF(i)/dble(NN)/RRR/RRR/coefficient
!end do
!close(26)
!
!
write(*,*) "ACP=",N_A,MC_max,dble(N_A)/dble(MC_max)
end program LJ_fluid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output(j,x)
use param
implicit none
integer,intent(in) :: j
double precision,dimension(dime,NN),intent(in) :: x
character(len=32) :: filename
integer :: i

write(filename,'("./data",i0,"/x_",i7.7,".dat")')ID,j
open(17,file=filename,status='replace')
do i=1,NN
   write(17,*) x(1,i),x(2,i),x(3,i)
end do
close(17)
end subroutine output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine random(IX,RND)
use param
implicit none
integer,intent(inout) :: IX
double precision,intent(out) :: RND
integer :: m,k 

!initial condition
data m,k/2147483647,48828125/  !m=2**31-1,k=5**11

!random number
IX=iand(IX*k,m)  !IX=IX*k ==> abs(IX)
RND=dble(IX)/dble(m)  !normalization

end subroutine random
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine energy(r,E_pair)
use param
implicit none
double precision,intent(in) :: r
double precision,intent(out) :: E_pair

E_pair=0.0d0
E_pair=4.0d0*(r**(-12.0d0)-r**(-6.0d0))
end subroutine energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine distance(x,n_i,n_j,r)
use param
implicit none
double precision,dimension(dime,NN),intent(in) :: x
integer,intent(in) :: n_i,n_j
double precision,intent(out) :: r
double precision :: rr
double precision,dimension(dime) :: x_ij
integer,dimension(dime) :: l
integer :: i,j

r=0.0d0
rr=0.0d0
x_ij=0.0d0
l=0.0d0

do i=1,dime
   x_ij(i)=x(i,n_i)-x(i,n_j)
   l(i)=int(2.0d0*x_ij(i)/length)
   x_ij(i)=x_ij(i)-l(i)*length

   rr=rr+x_ij(i)**2
end do

r=sqrt(rr)
end subroutine distance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initial(x)
use param
implicit none
double precision,dimension(dime,NN),intent(inout) :: x
integer :: i,j,xd

x=0.0d0


xd=int(NN**(1.0d0/3.0d0))+1

do i=1,NN
      x(1,i)=mod(i-1,xd)
      x(2,i)=mod((i-1)/xd,xd)
      x(3,i)=(i-1)/xd/xd

end do

x=length*x/dble(xd)

end subroutine initial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RDF_calc(x,step,CNT_RDF)
use param 
implicit none
double precision,dimension(dime,NN),intent(in) :: x
integer,intent(in) :: step
integer,dimension(N_hist),intent(inout) :: CNT_RDF
double precision :: r,rr,Ri,keisu
double precision,dimension(dime,NN*3**dime) :: x_add
double precision,dimension(dime,NN) :: x_new
double precision :: x_ij
integer,dimension(dime) :: l
integer,dimension(N_hist) :: CNT1
character(len=16) :: filenameRDF
integer :: i,j,k,ix,iy,iz,CNTa

do i=1,NN
   do j=1,dime
      l(j)=int(2*x(j,i)/length-1.0d0)
      x_new(j,i)=x(j,i)-l(j)*length
   end do
end do

do i=1,NN
   CNTa=1
   do ix=-1,1
   do iy=-1,1
   do iz=-1,1
      if(ix==0 .and. iy==0 .and. iz==0)then
         x_add(1,i)=x_new(1,i)+ix*length
         x_add(2,i)=x_new(2,i)+iy*length
         x_add(3,i)=x_new(3,i)+iz*length
      else 
         x_add(1,i+NN*CNTa)=x_new(1,i)+ix*length
         x_add(2,i+NN*CNTa)=x_new(2,i)+iy*length
         x_add(3,i+NN*CNTa)=x_new(3,i)+iz*length
         CNTa=CNTa+1
      end if
   end do
   end do
   end do
end do

CNT1=0   
do i=1,NN
   do j=1,NN*3**dime
      if(i/=j)then
         !!!!!!!!!distance calc!!!!!!!!!!
         rr=0.0d0
         do k=1,dime
            x_ij=x_add(k,i)-x_add(k,j)
            rr=rr+x_ij**2
         end do
         r=sqrt(rr)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         CNT1(int(r/dr)+1)=CNT1(int(r/dr)+1)+1
         CNT_RDF(int(r/dr)+1)=CNT_RDF(int(r/dr)+1)+1
      end if
   end do
end do

!write(filenameRDF,'("RDF_",i7.7,".dat")')step
!open(25,file="./RDF/"//filenameRDF,status='replace')
!keisu=4.0*pi*dble(NN)*dr/(length**3.0)
!do i=1,N_hist
!   Ri=dr*(dble(i)-1.0d0)
!   write(25,*) Ri,CNT_RDF(i)/dble(NN)/Ri/Ri/keisu
!end do
!close(25)

end subroutine RDF_calc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine outputDETAIL
use param
implicit none
character(len=32) file1


write(file1,'("./data",i0,"/Detail.dat")')ID
open(20,file=file1,status='replace')
write(20,*)"ID    :",ID
write(20,*)"seed  :",seed
write(20,*)"NN    :",NN
write(20,*)"MC_max:",MC_max
write(20,*)"eta   :",eta
write(20,*)"T     :",T
write(20,*)"random initial"
write(20,*)""

close(20)
end subroutine outputDETAIL
                                        
