function dim (x,y)
  
  use precision
  
  real(kind=realType) x,y,z
  real(kind=realType) :: dim

  dim = x - y 
  if (dim < 0.0) then
     dim = 0.0
  end if
         
end function dim
