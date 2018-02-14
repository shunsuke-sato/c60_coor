program main
  implicit none
  integer,parameter :: natom = 60
  real(8),parameter :: dist0 = 1d0
  real(8) :: rion(3,natom), rion0(3), rion_t(3)
  integer :: iatom,jatom,i
  real(8) :: r, diff_sum
  real(8) :: dist_r(3)
  integer :: dist_i(3)


  do iatom = 1,natom
    do i = 1,3
      call random_number(r)
      rion(i,iatom) = r
    end do
  end do

  rion0(1) = sum(rion(1,:))/natom
  rion0(2) = sum(rion(2,:))/natom
  rion0(3) = sum(rion(3,:))/natom

  rion(1,:) = rion(1,:) - rion0(1)
  rion(2,:) = rion(2,:) - rion0(2)
  rion(3,:) = rion(3,:) - rion0(3)


  do
    diff_sum = 0d0
    do iatom = 1,natom
      dist_r = 1d20
      dist_i = 0
      do jatom = 1,natom
        if(iatom == jatom)cycle
        r = sqrt( sum( (rion(:,iatom) - rion(:,jatom))**2 ))
        if(r < dist_r(1))then
          dist_r(3) = dist_r(2); dist_i(3) = dist_i(2)
          dist_r(2) = dist_r(1); dist_i(2) = dist_i(1)
          dist_r(1) = r        ; dist_i(1) = jatom
        else if(r < dist_r(2))then
          dist_r(3) = dist_r(2); dist_i(3) = dist_i(2)
          dist_r(2) = r        ; dist_i(2) = jatom
        else if(r < dist_r(3))then
          dist_r(3) = r        ; dist_i(3) = jatom
        end if

      end do

        if(abs(dist_r(1)-dist0) > max( abs(dist_r(2)-dist0), abs(dist_r(2)-dist0)))then
          i = 1
        else if (abs(dist_r(2)-dist0) > max( abs(dist_r(1)-dist0), abs(dist_r(3)-dist0)))then
          i = 2
        else 
          i = 3
        end if

        diff_sum = diff_sum + sum(abs(dist_r(1:3)-dist0))
       
        rion(:,iatom) = rion(:,iatom) &
          + dist0*(rion(:,dist_i(i)) - rion(:,iatom))/dist_r(dist_i(i))

    end do
    write(*,*)diff_sum
    if(diff_sum < 1d20)exit

  end do

  rion0(1) = sum(rion(1,:))/natom
  rion0(2) = sum(rion(2,:))/natom
  rion0(3) = sum(rion(3,:))/natom

  rion(1,:) = rion(1,:) - rion0(1)
  rion(2,:) = rion(2,:) - rion0(2)
  rion(3,:) = rion(3,:) - rion0(3)

  open(20,file="c60.xyz")
  do iatom = 1,natom
    write(20,"(I7,2x,999e26.16e3)")iatom,rion(:,iatom)
  end do
  close(20)


end program main

