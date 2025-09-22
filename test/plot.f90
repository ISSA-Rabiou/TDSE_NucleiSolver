program plot
   implicit none
 CHARACTER(len=100) :: name

 write(*,*) 'Quel est votre Nom svp ?'
 read(*,*) name
 write(*,*) 'Bonjour ', trim(name)
end program plot